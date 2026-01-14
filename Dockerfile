FROM python:3.11-slim

# Prevent Python from writing pyc files & buffering stdout
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# ---- System Dependencies ----
RUN apt-get update && apt-get install -y \
    procps \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# ---- python tooling ----
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Jupyter kernel support (THIS is what your IDE needs)
RUN pip install --no-cache-dir \
    ipykernel \
    jupyter-client \
    ipywidgets \
    tqdm

# ---- project workspace ----
WORKDIR /app

# Copy the source
COPY . .

# Register kernel explicitly (important for IDEs)
RUN python -m ipykernel install --sys-prefix --name deconomix --display-name "Python 3.11 (deconomix)"

# Default command: keep container alive for IDE attach
CMD ["sleep", "infinity"]
