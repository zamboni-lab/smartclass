FROM docker.io/library/python:3.12-slim
RUN apt update && apt install -y libxrender1 libxtst6 libxi6
RUN pip install uv
COPY uv.lock pyproject.toml ./
RUN uv sync
RUN mkdir /app
RUN adduser --system --no-create-home nonroot
USER nonroot
WORKDIR /app
# CMD [ "TODO"]
