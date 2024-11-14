## Run on Workstation

```bash
ssh -l pmurray pmurray.lin.flatironinstitute.org
cd ~/development/MutantAutoMate
git pull origin main
module load docker
docker build . -t mutantautomate
docker run -it --rm \
  -v $(pwd):/usr/src/app \
  -w /usr/src/app \
  -e FLASK_ENV=development \
  -p 5000:5000 \
  mutantautomate \
  flask --app src/mutantautomate/app2 run --host=0.0.0.0 --debug
docker ps
docker restart $(docker ps -q)
```