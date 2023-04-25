# crispr_microarray_mutant_significance_testing

## Preface

This was originally written as a programming test, though I thought the code would be better off published in case others wish to use this (especially between academia and students).  The original data I have permutated and applied random noise to, to yield a suitable example to execute this code with.


## Running Instructions

### Docker Image Building

Navigate to the folder of the docker file and run the following command, to generate the microarray_sig_test image

`docker build --tag microarray_sig_test . `


### Docker Running image to container

Execute the following command after building, to create and run a container called rest-server

`docker run --publish 8000:5000 --name rest-server microarray_sig_test `

Then execute this command to verify the output :

`curl localhost:8000`

alternatively, direct your favourite web browser to `localhost:8000`


For running other things (like the testing script), the bash shell can be opened:
`docker exec -t -i rest-server bash`  

running the testing script would be (within the /app directory)

`python -m unittest`

### Docker Termination

run this when finished with the container

`docker stop rest-server`

Likewise you may execute this to restart it (in the background)

`docker start rest-server`







## Methodology

This docker image was built to find whether the presence,as opposed to absence, of a mutation had an effect on the cell counts.  There are two data files - a table containing cell fold count data for each model (per gene knock-out), and a separate table telling us which mutations belonged to each cell model.  The script converts both tables to long format and merges the two tables, to yield a table of mutant-gene knockout pairs per model, whether the mutation was present, and their cell fold count change.

To assess the difference between presence/absence of a mutation for data per mutation-ko, The Brunner-Munzel non parametric test was used from scipy.  Whilst a more traditional choice would have been the Mann-Whitney-U test, the latter has more assumptions, notably that the sample variance is equivalent.  Thus, after calculating pairwise Brunner-Munzel tests, Benjamini-Hochberg False Discovery Rate p-value correction was applied, and selected those results with a corrected p value <0.05.  These results, would be where there is a “significant” difference between presence and absence.

### Criticisms

Significant doesn’t necessarily mean better.

![image](https://user-images.githubusercontent.com/37003581/234383650-3c82cf2c-7601-4571-8d2b-cde262e9760a.png)

In this histogram example with Gene4_mut mutation and GeneB_KO knockout, the orange bars show how the cell fold counts are distributed when the mutation is present, and the (very thin) blue bar shows where the absence of mutation data is distributed.  Arguably, the reason the distributions here show up as significant is due to poor sample size on behalf of the absence data (3 data points), thus realistically this is probably not a useful result and more data might yield a greater overlap.  
 
Thus, one may wish to cherry pick results to ensure that increasingly small sample sizes are avoided. 


