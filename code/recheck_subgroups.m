// Recheck: number of subgroups of (Z/N)^x containing -1, as enumerated by Subgroups(U) (the loop used by
// xdelta.m, xdelta_small.m, perdelta.m). Compared with Birkhoff's formula in code/recheck_subgroups.py.
SetColumns(0);
Ns := [13..500] cat [504, 512, 540, 567, 576, 600, 640, 648, 675, 720, 756, 768, 784, 800, 864, 900, 972, 1024, 1080, 1152, 1296, 1440, 1536, 1728, 2304, 2592, 3072, 3456, 4608, 6912];
for N in Ns do
  U, mU := UnitGroup(Integers(N)); m1 := (Integers(N)!(-1)) @@ mU;
  c := #[1 : rec in Subgroups(U) | m1 in rec`subgroup];
  printf "SUB N=%o count=%o\n", N, c;
end for;
print "DONE";
quit;
