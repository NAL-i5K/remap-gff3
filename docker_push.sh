echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin
docker tag remap-gff3 nali5k/remap-gff3
docker push nali5k/remap-gff3
