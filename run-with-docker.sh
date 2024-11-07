# docker build -t mutantautomate .
# module load docker
docker \
  run \
  --rm \
  --detach \
  --name mutantautomate \
  --volume $(pwd):/usr/src/app \
  --publish 5000:5000 \
  --workdir /usr/src/app \
  --env FLASK_ENV=development \
  --workdir /usr/src/app \
  mutantautomate \
  flask --app src/mutantautomate/app2 run --host=0.0.0.0 --debug