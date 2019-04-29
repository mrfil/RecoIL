!mex -v LDOPTIMFLAGS="-O2" CXXOPTIMFLAGS="-O2 -DNDEBUG -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" interp2_table_mex.cpp interp2_table1_for.cpp
!mex -v LDOPTIMFLAGS="-O2" CXXOPTIMFLAGS="-O2 -DNDEBUG -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" interp2_table_adj_mex.cpp interp2_table1_adj.cpp

%!mex -v CFLAGS="\$CFLAGS " LDFLAGS="\$LDFLAGS " interp2_table_mex.cpp interp2_table1_for.cpp
%!mex -v CFLAGS="\$CFLAGS " LDFLAGS="\$LDFLAGS " interp2_table_adj_mex.cpp interp2_table1_adj.cpp