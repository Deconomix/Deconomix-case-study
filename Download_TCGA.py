#!/usr/bin/env python3
import os, re, io, json, gzip, tarfile, argparse
from typing import List, Optional, Tuple
import requests
import pandas as pd
import numpy as np
from scipy import sparse
import openpyxl
from tqdm import tqdm

try:
    import anndata as ad
    HAVE_ANNDATA = True
except Exception:
    HAVE_ANNDATA = False

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"


from typing import Optional

def download_brca_pam50_metadata() -> Optional[pd.DataFrame]:
    """
    Download BRCA PAM50 metadata from the supplementary Excel file.
    
    Returns:
    --------
    pd.DataFrame or None
        DataFrame with Sample_ID and BRCA_Subtype_PAM50 columns
    """
    
    url = "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301193-mmc4.xlsx"
    
    try:
        print("Downloading BRCA PAM50 metadata...")
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        print(f"Downloaded {len(response.content)} bytes")
        
        # Read the Excel file
        excel_data = pd.read_excel(io.BytesIO(response.content), sheet_name=None, skiprows=1)
        
        print(f"Excel file contains {len(excel_data)} sheets:")
        for sheet_name in excel_data.keys():
            print(f"  - {sheet_name}")
        
        # Try to find the sheet with PAM50 data
        pam50_df = None
        
        for sheet_name, df in excel_data.items():
            print(f"\nChecking sheet '{sheet_name}':")
            print(f"  Shape: {df.shape}")
            print(f"  Columns: {list(df.columns)[:10]}...")  # Show first 10 columns
            
            # Look for Sample_ID and BRCA_Subtype_PAM50 columns
            sample_id_cols = [col for col in df.columns if 'sample' in str(col).lower() and 'id' in str(col).lower()]
            pam50_cols = [col for col in df.columns if 'brca_subtype_pam50' in str(col).lower()]
            
            if sample_id_cols and pam50_cols:
                print(f"  Found Sample_ID columns: {sample_id_cols}")
                print(f"  Found PAM50 columns: {pam50_cols}")
                
                # Extract the relevant columns
                sample_col = sample_id_cols[0]
                pam50_col = pam50_cols[0]
                
                pam50_df = df[[sample_col, pam50_col]].copy()
                pam50_df.columns = ['Sample_ID', 'BRCA_Subtype_PAM50']
                
                # Remove any rows with missing values
                pam50_df = pam50_df.dropna()
                
                print(f"  Extracted {len(pam50_df)} samples with PAM50 data")
                break
            
            # Alternative: look for any column containing 'pam50'
            elif any('pam50' in str(col).lower() for col in df.columns):
                pam50_related_cols = [col for col in df.columns if 'pam50' in str(col).lower()]
                print(f"  Found PAM50-related columns: {pam50_related_cols}")
                
                # Try to find a sample ID column
                potential_id_cols = [col for col in df.columns if 'sample' in str(col).lower() or 'id' in str(col).lower() or 'tcga' in str(col).lower()]
                
                if potential_id_cols:
                    print(f"  Potential ID columns: {potential_id_cols}")
                    
                    # Use the first PAM50 column and first ID column
                    id_col = potential_id_cols[0]
                    pam50_col = pam50_related_cols[0]
                    
                    pam50_df = df[[id_col, pam50_col]].copy()
                    pam50_df.columns = ['Sample_ID', 'BRCA_Subtype_PAM50']
                    pam50_df = pam50_df.dropna()
                    
                    print(f"  Extracted {len(pam50_df)} samples using alternative method")
                    break
        
        if pam50_df is not None:
            print(f"\nSuccessfully extracted PAM50 data:")
            print(f"  Total samples: {len(pam50_df)}")
            print(f"  Unique subtypes: {pam50_df['BRCA_Subtype_PAM50'].nunique()}")
            print(f"  Subtype distribution:")
            print(pam50_df['BRCA_Subtype_PAM50'].value_counts().to_string())
            
            return pam50_df
        else:
            print("\nCould not find Sample_ID and BRCA_Subtype_PAM50 columns in any sheet")
            
            # Show all available columns for debugging
            print("\nAll columns across all sheets:")
            for sheet_name, df in excel_data.items():
                print(f"\nSheet '{sheet_name}' columns:")
                for i, col in enumerate(df.columns):
                    print(f"  {i+1:2d}. {col}")
            
            return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        return None
    except Exception as e:
        print(f"Error processing Excel file: {e}")
        return None

def save_pam50_data(df: pd.DataFrame, filename: str = "brca_pam50_subtypes.csv") -> None:
    """
    Save the PAM50 data to a CSV file.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with Sample_ID and BRCA_Subtype_PAM50
    filename : str
        Output filename
    """
    try:
        df.to_csv(filename, index=False)
        print(f"PAM50 data saved to: {filename}")
    except Exception as e:
        print(f"Error saving file: {e}")

def convert_sample_id_to_patient_id(sample_id: str) -> str:
    """
    Convert TCGA sample ID to patient ID format.
    
    TCGA sample IDs are typically: TCGA-XX-XXXX-XXX
    Patient IDs are typically: TCGA-XX-XXXX
    
    Parameters:
    -----------
    sample_id : str
        TCGA sample ID
        
    Returns:
    --------
    str
        TCGA patient ID
    """
    if pd.isna(sample_id):
        return sample_id
    
    # Split by hyphen and take first 3 parts
    parts = str(sample_id).split('-')
    if len(parts) >= 3:
        return '-'.join(parts[:3])
    else:
        return str(sample_id)

def process_pam50_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process the PAM50 data to add patient IDs and clean up.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Raw PAM50 data with Sample_ID and BRCA_Subtype_PAM50
        
    Returns:
    --------
    pd.DataFrame
        Processed DataFrame with additional patient ID column
    """
    df_processed = df.copy()
    
    # Add patient ID column
    df_processed['Patient_ID'] = df_processed['Sample_ID'].apply(convert_sample_id_to_patient_id)
    
    # Clean up sample IDs (remove any whitespace)
    df_processed['Sample_ID'] = df_processed['Sample_ID'].str.strip()
    df_processed['Patient_ID'] = df_processed['Patient_ID'].str.strip()
    
    # Clean up subtype names
    df_processed['BRCA_Subtype_PAM50'] = df_processed['BRCA_Subtype_PAM50'].str.strip()
    
    # Reorder columns
    df_processed = df_processed[['Sample_ID', 'Patient_ID', 'BRCA_Subtype_PAM50']]
    
    print(f"\nProcessed PAM50 data:")
    print(f"  Columns: {list(df_processed.columns)}")
    print(f"  Sample of data:")
    print(df_processed.head().to_string())
    
    return df_processed


def _detect_workflow(project_id: str) -> str:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["RNA-Seq"]}},
        ],
    }
    params = {"filters": json.dumps(filters), "size": 0, "facets": "analysis.workflow_type"}
    r = requests.get(GDC_FILES, params=params, timeout=120)
    r.raise_for_status()
    agg = r.json()["data"].get("aggregations", {})
    buckets = (agg.get("analysis.workflow_type") or {}).get("buckets", [])
    if not buckets:
        raise RuntimeError("No RNA-Seq Gene Expression Quantification files found.")
    available = [b["key"] for b in buckets]
    # preference order
    for w in ["STAR - Counts", "HTSeq - Counts", "HTSeq - FPKM-UQ", "HTSeq - FPKM"]:
        if w in available:
            return w
    return available[0]

def _gdc_query_files(project_id: str, workflow_type: Optional[str] = None, size: int = 2000) -> pd.DataFrame:
    if workflow_type is None:
        workflow_type = _detect_workflow(project_id)
        print(f"Using workflow_type: {workflow_type}")

    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["RNA-Seq"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": [workflow_type]}},
        ],
    }
    fields = [
        "file_id","file_name","md5sum","file_size","state","data_format","created_datetime",
        "analysis.workflow_type","data_category","data_type","experimental_strategy",
        "cases.project.project_id","cases.disease_type","cases.submitter_id",
        "cases.samples.submitter_id","cases.samples.sample_type",
        "cases.samples.portions.analytes.aliquots.submitter_id",
        "associated_entities.entity_type","associated_entities.entity_submitter_id",
    ]
    params = {"filters": json.dumps(filters), "fields": ",".join(fields), "format": "JSON", "size": size}
    r = requests.get(GDC_FILES, params=params, timeout=180)
    r.raise_for_status()
    hits = r.json()["data"].get("hits", [])

    if not hits:
        raise RuntimeError("No files returned by GDC for the chosen filters.")
    rows = []
    for h in hits:
        file_uuid = h.get("file_id")
        cases = h.get("cases") or []
        case_submitter_id = None
        case_project = None
        disease_type = None
        sample_submitter_id = None
        sample_type = None
        aliquot_from_assoc = None
        for ae in h.get("associated_entities") or []:
            if ae.get("entity_type") == "aliquot":
                aliquot_from_assoc = ae.get("entity_submitter_id")
                break
        if cases:
            c0 = cases[0]
            case_submitter_id = c0.get("submitter_id")
            proj = c0.get("project") or {}
            case_project = proj.get("project_id")
            disease_type = c0.get("disease_type")
            samples = c0.get("samples") or []
            if aliquot_from_assoc:
                found = False
                for s in samples:
                    aliquots = []
                    for portion in s.get("portions") or []:
                        for analy in portion.get("analytes") or []:
                            for aq in analy.get("aliquots") or []:
                                sid = aq.get("submitter_id")
                                if sid: aliquots.append(sid)
                    if aliquot_from_assoc in aliquots:
                        sample_submitter_id = s.get("submitter_id")
                        sample_type = s.get("sample_type")
                        found = True
                        break
                if not found and samples:
                    s0 = samples[0]; sample_submitter_id = s0.get("submitter_id"); sample_type = s0.get("sample_type")
            elif samples:
                s0 = samples[0]; sample_submitter_id = s0.get("submitter_id"); sample_type = s0.get("sample_type")
        rows.append(dict(
            file_id=file_uuid,
            file_name=h.get("file_name"),
            md5sum=h.get("md5sum"),
            file_size=h.get("file_size"),
            state=h.get("state"),
            data_format=h.get("data_format"),
            created_datetime=h.get("created_datetime"),
            workflow_type=(h.get("analysis") or {}).get("workflow_type"),
            data_category=h.get("data_category"),
            data_type=h.get("data_type"),
            experimental_strategy=h.get("experimental_strategy"),
            project_id=case_project,
            disease_type=disease_type,
            case_submitter_id=case_submitter_id,
            sample_submitter_id=sample_submitter_id,
            sample_type=sample_type,
            aliquot_submitter_id=aliquot_from_assoc,
        ))
    df = pd.DataFrame(rows)
  
    return df


def get_tcga_rnaseq(
    project_id: str = "TCGA-BRCA",
    workflow_type: Optional[str] = None,
    out_dir: str = "./gdc_brca",
):
    os.makedirs(out_dir, exist_ok=True)
    df = _gdc_query_files(project_id=project_id, workflow_type=workflow_type)
    
    file_ids = df["file_id"].tolist()
    file_ids = [str(fid) for fid in file_ids]

    # Check if out_dir ("gdc_brca" by default) is empty
    if os.listdir(out_dir):
        print(f"Directory '{out_dir}' is not empty. Skipping download.")
    else:
        for i in range(0, len(file_ids), 100):
            print(f"Downloading files {i+1} to {min(i+100, len(file_ids))} of {len(file_ids)}")
            batch = file_ids[i:i+100]

            params = {"ids": batch}

            response = requests.post(GDC_DATA,
                                data = json.dumps(params),
                                headers={
                                    "Content-Type": "application/json"
                                    })

            response_head_cd = response.headers["Content-Disposition"]

            file_name = re.findall("filename=(.+)", response_head_cd)[0]

            with open(os.path.join(out_dir, file_name), "wb") as output_file:
                output_file.write(response.content)

    # Downloading PAM50 metadata
    pam50_df = download_brca_pam50_metadata()
    pam50_df.set_index("Sample_ID", inplace=True)

    if pam50_df is not None:
        
        # Show some statistics
        print(f"\nFinal statistics:")
        print(f"  Total samples: {len(pam50_df)}")
        print(f"  Unique patients: {pam50_df.index.nunique()}")
        print(f"  PAM50 subtype distribution:")
        subtype_counts = pam50_df['BRCA_Subtype_PAM50'].value_counts()
        for subtype, count in subtype_counts.items():
            percentage = (count / len(pam50_df)) * 100
            print(f"    {subtype}: {count} ({percentage:.1f}%)")
            
    tar_files = [os.path.join(out_dir, f) for f in os.listdir(out_dir) if f.endswith(".tar.gz")]

    all_data = []
    pam50_list = []

    for tf in tar_files:
        print(f"Extracting {tf}...")
        tf_dir = tf[:-7]
        # _extract_tar(tf, tf_dir)

        manifest_path = os.path.join(tf_dir, "MANIFEST.txt")
        manifest = pd.read_csv(manifest_path, sep="\t")
        ids = manifest["id"].tolist()
        file_names = manifest["filename"].tolist()
        for f, fid in zip(file_names, ids):
            case_id = None
            sample_type = None
            for idx, row in df.iterrows():
                if row["file_id"] == fid:
                    case_id = row["case_submitter_id"]
                    sample_type = row["sample_type"]
                    break
            PAM50 = pam50_df.loc[case_id]["BRCA_Subtype_PAM50"] if (pam50_df is not None and case_id in pam50_df.index) else "Unknown"
            pam50_list.append(PAM50)

            print(f"File: {f}, File ID: {fid}, Case Submitter ID: {case_id}")
            counts = pd.read_csv(os.path.join(tf_dir, f), sep="\t", skiprows=1)
            counts = counts.iloc[4:,:]


            counts_df = counts.set_index("gene_name")[["unstranded"]].copy()
            counts_df.columns = [case_id]
            counts_df.loc["sample_type"] = [sample_type]
            counts_df = counts_df.loc[["sample_type"] + [idx for idx in counts_df.index if idx != "sample_type"]]
            counts_df = counts_df[~counts_df.index.duplicated(keep='first')]
            all_data.append(counts_df)

    
    combined_df = pd.concat(all_data, axis=1, join='inner')
    print(f"Nan values in combined data: {combined_df.isna().sum().sum()}")

    combined_df.loc["PAM50"] = pam50_list
    combined_df = combined_df.loc[["PAM50"] + [idx for idx in combined_df.index if idx != "sample_type" and idx != "PAM50"]]

    combined_df.to_csv("tcga_brca_rnaseq_counts.csv")



def main():
    get_tcga_rnaseq(
        project_id="TCGA-BRCA",
        workflow_type=None,
        out_dir="./Data/TCGA",
    )

    print("Done.")

if __name__ == "__main__":
    main()