# Deployment

## Pre-requisites

* sudo apt-get install docker docker-compose

* Download a linux yasara installation executable
  (see: http://www.yasara.org/update.htm)
  to the root directory and name it install_yasara_16.4.6

## Installation

1. In docker-compose.yml, make sure that the volume /mnt/cmbi4 mounts to a
   location where dssp files are stored. At the same time, make the volume
   /data point to a location with sufficient disk space to store models and
   blast databanks.
2. Run 'docker-compose up' from the root directory or for testing purposes:
   'docker-compose -f docker-compose.yml -f docker-compose-dev.yml'
