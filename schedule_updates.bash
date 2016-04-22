#!/bin/bash

MYDIR=$(dirname $0)

UPDATE=$MYDIR/update_databanks.bash

chmod 755 $UPDATE

if ! [ -f /data/blast/uniprot.pal ] || ! [ -f /data/blast/templates.psq ] || ! [-f /data/blacklisted_templates ]; then

    # Databanks not present, build now!
    /bin/bash $UPDATE
fi

# Copy this script's environment variables to cron job:
ENV="
PATH=$PATH
"

CRONFILE=$MYDIR/update_cron

/bin/echo -e "$ENV\n0 20 * * 5 /bin/bash $UPDATE\n" > $CRONFILE
crontab $CRONFILE

cron -f
