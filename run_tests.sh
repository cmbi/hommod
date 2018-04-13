#!/usr/bin/env bash

while [[ $# -gt 1 ]]
do
_key="$1"

case $_key in
    -t|--tests)
        _tests="$2"
        shift
        ;;
    *)
        echo "Unknown option: $_key"
        exit 1
        ;;
esac
shift # past argument or value
done

if [ "x$_tests" = "x" ]; then
    _tests=tests/unit
fi

echo $_tests

case "$_tests" in
    *unit*)
        _dc_run_opts="--no-deps --rm"
        ;;
    *)
        _dc_run_opts="--rm"
        ;;
esac

_nose_opts="--with-coverage --cover-inclusive --cover-package hommod"
_dc_opts="-f docker-compose.yml -f docker-compose-dev.yml"
_command="docker-compose $_dc_opts run $_dc_run_opts celery nosetests $_nose_opts $_tests -v"
echo $_command
$_command
exit_code=$?

# Remove all containers and network only if the tests passed. Keep them around
# for debugging the failed tests.
if [ $exit_code -eq 0 ]; then
    docker-compose down
fi
exit $exit_code
