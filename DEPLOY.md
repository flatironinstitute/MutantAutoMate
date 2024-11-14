## Run on Workstation

```bash
ssh -l pmurray pmurray.lin.flatironinstitute.org # login to workstation
cd ~/development/MutantAutoMate
git pull origin main
module load docker
docker build . -t mutantautomate
bash run-with-docker.sh
docker ps
docker restart $(docker ps -q)
```