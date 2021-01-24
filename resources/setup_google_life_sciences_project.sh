# This script is sourced from
# - https://cloud.google.com/life-sciences/docs/tutorials/nextflow
# - https://www.nextflow.io/docs/latest/google.html

# Setup gcloud CLI


export PROJECT=bahia-analysis

gcloud projects create $PROJECT --name=$PROJECT

gcloud config set project bahia-analysis

export BUCKET=bahia-analysis-bucket

gsutil mb gs://$BUCKET

export SERVICE_ACCOUNT_NAME=nextflow-service-account
export SERVICE_ACCOUNT_ADDRESS=${SERVICE_ACCOUNT_NAME}@${PROJECT}.iam.gserviceaccount.com


gcloud iam service-accounts create ${SERVICE_ACCOUNT_NAME}

gcloud projects add-iam-policy-binding ${PROJECT} \
    --member serviceAccount:${SERVICE_ACCOUNT_ADDRESS} \
    --role roles/lifesciences.workflowsRunner

gcloud projects add-iam-policy-binding ${PROJECT} \
    --member serviceAccount:${SERVICE_ACCOUNT_ADDRESS} \
    --role roles/iam.serviceAccountUser

gcloud projects add-iam-policy-binding ${PROJECT} \
    --member serviceAccount:${SERVICE_ACCOUNT_ADDRESS} \
    --role roles/serviceusage.serviceUsageConsumer

gcloud projects add-iam-policy-binding ${PROJECT} \
    --member serviceAccount:${SERVICE_ACCOUNT_ADDRESS} \
    --role roles/storage.objectAdmin


export SERVICE_ACCOUNT_KEY=${SERVICE_ACCOUNT_NAME}-private-key.json

gcloud iam service-accounts keys create \
  --iam-account=${SERVICE_ACCOUNT_ADDRESS} \
  --key-file-type=json ${SERVICE_ACCOUNT_KEY}

export SERVICE_ACCOUNT_KEY_FILE=${PWD}/${SERVICE_ACCOUNT_KEY}

export GOOGLE_APPLICATION_CREDENTIALS=${PWD}/${SERVICE_ACCOUNT_KEY}

nextflow run main.nf -params-file params.yml -profile google_cloud_life_sciences -queue-size 1



