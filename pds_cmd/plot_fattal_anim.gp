plot sprintf("pnt_%d_%d.tsv", q,i) using (2**(q+1))*($1):(2**(q+1))*($2)
pause -1
i=i+1
if(i<30) reread
