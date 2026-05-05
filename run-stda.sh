# Create a directory for xyz files
XYZ_DIR=crystal_geom_xyz
mkdir -p $XYZ_DIR
# Use openbabel to convert mol files into xyz files
MOL_DIR=crystal_geom
for INFILE in $MOL_DIR/*.mol
do
    BASENAME=$(basename $INFILE)
    MOLID=${BASENAME%.mol}
    OUTFILE=$XYZ_DIR/$MOLID.xyz
    obabel -imol $INFILE -oxyz -O$OUTFILE
done

LOG_DIR=crystal_geom_stda_logs
# Move to the directory containing the optimized xyz files
# Run xtb and stda on each one, and save the log file containing the energy
mkdir -p $LOG_DIR
cd $XYZ_DIR
for INFILE in *.xyz
do
    MOLID=${INFILE%.xyz}
    xtb4stda $INFILE
    stda -xtb > ../$LOG_DIR/$MOLID.log
done
cd ..

uv run summarize_stda.py $LOG_DIR crystal_stda.csv.gz
Rscript compare_stda.R
