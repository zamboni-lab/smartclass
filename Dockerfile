FROM python:3.12.0-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxi6 \
    libxrender1 \
    libxtst6 \
    && rm -rf /var/lib/apt/lists/*
RUN pip install uv
COPY uv.lock pyproject.toml ./
RUN --mount=type=cache,target=/root/.cache/pip uv sync
RUN mkdir /app
RUN adduser --system --no-create-home nonroot
USER nonroot
WORKDIR /app
# CMD [ "TODO"]
