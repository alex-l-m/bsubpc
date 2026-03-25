mkdir -p 3d_diagrams

for INPATH in crystal_geom/*.mol
do
    MOLID=`basename $INPATH .mol`
    python render_mol.py $INPATH 3d_diagrams/${MOLID}.png
done
