# docker build -t mutantautomate .
module load docker
docker run -it --rm \
  -v $(pwd):/usr/src/app \
  -w /usr/src/app \
  -e FLASK_ENV=development \
  -p 5000:5000 \
  mutantautomate \
  flask --app src/mutantautomate/app2 run --host=0.0.0.0 --debug