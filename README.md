# `mtb_mono_inh` nextflow pipeline
A pipeline for analysing _Mycobacterium tuberculosis_ samples.

## Minimal requirements (for local execution)

* Nextflow VERSION > 21.04
* Java >= 8
* Docker 

## Pipeline workflow

![dag file](./resources/dag.png)

## Quick start

### Local execution
1. Install nextflow 

	Please refer to [Nextflow page on github](https://github.com/nextflow-io/nextflow/) for more info.

2. Run it!

```
	nextflow run https://github.com/bioinformatics-lab/mtb_mono_inh_nf -params-file YOUR_PARAMS_FILE

```

$YOUR_PARAMS_FILE = file path, please refer to params/params.yml as a template to create your own params file

## Configuration Profiles.

You can use diferent profiles for this pipeline, based on the computation enviroment at your disposal contained within `conf`folder.

**Note**: Update conf/profile with your own configs.

## Tower execution
This Pipeline can be launched on `Tower`, please refer to [Tower launch documentation](https://help.tower.nf) for step-by-step execution tutorial.

When launching from `Tower`, please update and use the `params.yml` file contents.

## Mock execution using stub-run
This project has the `-stub-run` feature, that can be used for testing propouse, it can be used on `Tower` with the Advanced settings on launch. You can also test it locally, using the following command:

```
cd data/mock_data/ && bash generate_mock_files.sh

nextflow run main.nf \
		 -profile conf/stub.config \
		 -stub-run
``` 
