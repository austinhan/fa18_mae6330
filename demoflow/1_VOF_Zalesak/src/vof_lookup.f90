module vof_lookup
  use demoflow
  implicit none

  ! Number of new vertices on cut plane
  integer, dimension(16) :: cut_nvert
  data cut_nvert(1:16) / 0, 3, 3, 4, 3, 4, 4, 3, 3, 4, 4, 3, 4, 3, 3, 0/

  ! Number of resulting tets
  integer, dimension(16) :: cut_ntets
  data cut_ntets(1:16) / 1, 4, 4, 6, 4, 6, 6, 4, 4, 6, 6, 4, 6, 4, 4, 1/

  ! Index of first positive tet = # tets - # negative tets + 1
  integer, dimension(16) :: cut_nntet
  data cut_nntet(1:16) / 1, 2, 2, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2/
  
  ! First point on intersection 
  integer, dimension(4,16) :: cut_v1
  data cut_v1(1:4, 1) /-1,-1,-1,-1/
  data cut_v1(1:4, 2) / 1, 1, 1,-1/
  data cut_v1(1:4, 3) / 2, 2, 2,-1/
  data cut_v1(1:4, 4) / 1, 2, 1, 2/
  data cut_v1(1:4, 5) / 3, 3, 3,-1/
  data cut_v1(1:4, 6) / 1, 3, 1, 3/
  data cut_v1(1:4, 7) / 2, 3, 2, 3/
  data cut_v1(1:4, 8) / 4, 4, 4,-1/
  data cut_v1(1:4, 9) / 4, 4, 4,-1/
  data cut_v1(1:4,10) / 1, 4, 1, 4/
  data cut_v1(1:4,11) / 2, 4, 2, 4/
  data cut_v1(1:4,12) / 3, 3, 3,-1/
  data cut_v1(1:4,13) / 3, 4, 3, 4/
  data cut_v1(1:4,14) / 2, 2, 2,-1/
  data cut_v1(1:4,15) / 1, 1, 1,-1/
  data cut_v1(1:4,16) /-1,-1,-1,-1/

  ! Second point on intersection 
  integer, dimension(4,16) :: cut_v2
  data cut_v2(1:4, 1) /-1,-1,-1,-1/
  data cut_v2(1:4, 2) / 2, 3, 4,-1/
  data cut_v2(1:4, 3) / 3, 4, 1,-1/
  data cut_v2(1:4, 4) / 4, 4, 3, 3/
  data cut_v2(1:4, 5) / 4, 1, 2,-1/
  data cut_v2(1:4, 6) / 4, 4, 2, 2/
  data cut_v2(1:4, 7) / 4, 4, 1, 1/
  data cut_v2(1:4, 8) / 1, 2, 3,-1/
  data cut_v2(1:4, 9) / 1, 2, 3,-1/
  data cut_v2(1:4,10) / 3, 3, 2, 2/
  data cut_v2(1:4,11) / 3, 3, 1, 1/
  data cut_v2(1:4,12) / 4, 1, 2,-1/
  data cut_v2(1:4,13) / 2, 2, 1, 1/
  data cut_v2(1:4,14) / 3, 4, 1,-1/
  data cut_v2(1:4,15) / 2, 3, 4,-1/
  data cut_v2(1:4,16) /-1,-1,-1,-1/
  
  ! Vertices in each tet
  integer, dimension(4,6,16) :: cut_vtet
  data cut_vtet(1:4,1:6, 1) / 1, 2, 3, 4, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 2) / 5, 7, 6, 1,  6, 2, 3, 4,  4, 2, 5, 6,  5, 6, 7, 4, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 3) / 7, 5, 6, 2,  1, 3, 4, 6,  1, 5, 3, 6,  5, 7, 6, 1, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 4) / 5, 8, 6, 2,  5, 7, 8, 1,  5, 1, 8, 2,  5, 6, 8, 4,  5, 8, 7, 3,  5, 8, 3, 4/
  data cut_vtet(1:4,1:6, 5) / 6, 5, 7, 3,  2, 1, 4, 6,  6, 5, 4, 2,  6, 7, 5, 2, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 6) / 5, 6, 8, 3,  5, 8, 7, 1,  5, 8, 1, 3,  5, 8, 6, 4,  5, 7, 8, 2,  5, 8, 4, 2/
  data cut_vtet(1:4,1:6, 7) / 8, 6, 5, 3,  5, 7, 8, 2,  8, 5, 2, 3,  8, 5, 6, 4,  5, 8, 7, 1,  5, 8, 1, 4/
  data cut_vtet(1:4,1:6, 8) / 1, 2, 3, 7,  1, 2, 7, 6,  5, 7, 6, 1,  5, 6, 7, 4, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 9) / 5, 6, 7, 4,  1, 2, 3, 6,  5, 1, 3, 6,  5, 7, 6, 3, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,10) / 5, 8, 6, 4,  5, 7, 8, 1,  5, 8, 4, 1,  5, 6, 8, 3,  5, 8, 7, 2,  5, 8, 2, 3/
  data cut_vtet(1:4,1:6,11) / 8, 5, 6, 4,  5, 8, 7, 2,  8, 2, 5, 4,  8, 6, 5, 3,  5, 7, 8, 1,  5, 8, 3, 1/
  data cut_vtet(1:4,1:6,12) / 1, 4, 2, 7,  4, 1, 6, 7,  6, 7, 5, 4,  6, 5, 7, 3, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,13) / 8, 6, 5, 4,  5, 7, 8, 3,  8, 4, 5, 3,  8, 5, 6, 2,  5, 8, 7, 1,  5, 8, 1, 2/
  data cut_vtet(1:4,1:6,14) / 3, 4, 1, 7,  7, 6, 3, 4,  7, 6, 5, 3,  7, 5, 6, 2, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,15) / 7, 4, 2, 3,  2, 3, 6, 7,  5, 6, 7, 2,  5, 7, 6, 1, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,16) / 1, 2, 3, 4, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1/

  ! Side of cut plane (used to update i,j,k)
  integer, dimension(6,  16) :: cut_side
  data cut_side(1:6, 1) / 1,-1,-1,-1,-1,-1/
  data cut_side(1:6, 2) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6, 3) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6, 4) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6, 5) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6, 6) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6, 7) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6, 8) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6, 9) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6,10) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6,11) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6,12) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6,13) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6,14) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6,15) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6,16) / 2,-1,-1,-1,-1,-1/

  ! Vertices on cell to 6 tets for pure semi-Lagrangian
  integer, dimension(4,6) :: tet_map2
  data tet_map2(1:4,1) / 7, 4, 3, 6 /
  data tet_map2(1:4,2) / 6, 3, 2, 4 /
  data tet_map2(1:4,3) / 6, 2, 1, 4 /
  data tet_map2(1:4,4) / 7, 8, 4, 6 /
  data tet_map2(1:4,5) / 6, 5, 8, 4 /
  data tet_map2(1:4,6) / 6, 5, 4, 1 /
  
end module vof_lookup
