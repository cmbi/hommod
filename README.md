# Deployment

## Pre-requisites

* sudo apt-get install docker docker-compose

* Download a linux yasara installation executable
  (see: http://www.yasara.org/update.htm)
  to the root directory and name it install_yasara_<version>

## Installation

1. In docker-compose.yml, make sure that the volume /mnt/chelonium mounts to a
   location where dssp files are stored. At the same time, make the volume
   /data point to a location with sufficient disk space to store models and
   blast databanks.
2. Build the images: 'docker-compose build'
3. The databanks update script will run periodically. However if you start first time,
   you must build the databanks using: 'docker-compose run databanks ./update_databanks.bash'
4. Run 'docker-compose up' from the root directory or for testing purposes:
   ./run_dev.sh
