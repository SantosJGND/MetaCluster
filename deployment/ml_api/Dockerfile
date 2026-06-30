FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN apt-get update && \
    apt-get install -y git && \
    rm -rf /var/lib/apt/lists/*


ENV MLFLOW_TRACKING_URI=http://mlflow:5000
ENV MLFLOW_EXPERIMENT=INSaFLU_ML_API

EXPOSE 8000

CMD ["uvicorn", "ml_api.app:app", "--host", "0.0.0.0", "--port", "8000"]
