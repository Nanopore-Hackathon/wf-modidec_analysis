# wf-modidec (Part3: Analysis)
An Epi2Me integrated tool to detect and classify RNA modifications. 


## Introduction
In this directory you find part 3 of the code from ModiDec. The scripts located in ./bin facilitate the application of the trained model. You can find test data in the ./test folder of this repository.


### Functionality Overview
Below is a graphical overview of suggested routes through the pipeline depending on the desired output.

[image]

## Quickstart
Get your system up and ready
    - Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.10.0`)
    - Install [`Miniconda`](https://conda.io/miniconda.html)
    - Install [`Docker`](https://conda.io/miniconda.html)
    - Install [`Epi2Me Desktop`](https://labs.epi2me.io) (v5.1.14 or later)
    - Clone the Github repository (we recommend the GitHub Desktop client)
    - Clone ONT's kmer tabel github repository top obtain necessary [`kmer-tables`](https://github.com/nanoporetech/kmer_models) 

Import the workflow in Epi2Me without having a public repository
   -  Download the repository as a zip file
   -  Navigate in your filesystem to ~/epi2melabs/workflows
   -  Make a "modidec" directory
   -  Unpack your repository in ~/epi2melabs/workflows/modidec

Open Epi2Me
   -  Continue without signing in
   -  Click on View workflows
   -  You should now see wf-modidec_analysis appear as a workflow (Here we often experience some trouble, so feel free to inform us if something is wrong)
   -  Try to start the workflow and explore the menu of the workflow (The workflow will probably crash due to the lack of GPUs in your device) 
   -  Please have another look into the nextflow_schema.json file for a better understanding of how the menu is built.



## Guideline
Key aim: This script is already compatible with Epi2Me and can be explored as a show case. 


1. Navigate throug the functions and logic of Nextflow processes and the Nextflow pipeline of this repository. Try to understand the logics of the three scripts main.nf, nextflow.config and nextflow_schema.json
    
2. Check out the config.yaml file to see which input parameters have been necessary to run the pipeline.

3. Open the bin folder and explore the analysis_neural_network.py file
   -  Check out how the argument parser was defined and how the python file is structured.

4. Navigate to the envs folder and investigate the enviroment.yaml and the Dockerfile
   -  Find out how the environment.yaml is used to build the docker container
   -  Try to build the docker container yourself
   -  You can also try to pull our prebuilt docker container "stegiopast/modidec" from docker hub
   -  The container can be used for all three parts of this hackathon


## Credits & License

This code is provided by Dr.Nicolo Alagna and the Computational Systems Genetics Group of the University Medical Center of Mainz. © 2024 All rights reserved.

For the purpose of the Nanopore Hackathon, participants and collaborators of the event are granted a limited, non-exclusive, non-transferable license to use and modify the Applications for their intended purpose during the Time of the Event. Any further unauthorized reproduction, distribution, modification or publication of the Applications or their contents violates intellectual property rights and is strictly prohibited. For more please see the [Terms and Conditions](https://drive.google.com/file/d/18WN3YRoY9YvpYq6RCtwUQre-VAbN7jH6/view?usp=sharing) of the event.


