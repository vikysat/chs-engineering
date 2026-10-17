load data/pdb1cgk.ent, ref
load data/models/Q8H4L3_aligned.pdb, rice
load data/models/P48398_aligned.pdb, ipo
hide everything
bg_color white
set ray_opaque_background, on

select refA, ref and chain A and polymer
show cartoon, refA
color grey80, refA
show cartoon, rice
color cyan, rice
show cartoon, ipo
color orange, ipo
set cartoon_transparency, 0.5

select cys, refA and resi 164 and name SG
select pk_ref,  byres (refA and polymer) within 8 of cys
select pk_rice, byres rice within 8 of cys
select pk_ipo,  byres ipo within 8 of cys
show sticks, pk_ref or pk_rice or pk_ipo
color grey50, pk_ref and elem C
color cyan,   pk_rice and elem C
color orange, pk_ipo and elem C
show spheres, cys
set sphere_scale, 0.5
color yellow, cys

zoom pk_ref, 6
set stick_radius, 0.15
set ray_shadows, off
ray 1600, 1200
png figures/pocket_overlay.png
quit
