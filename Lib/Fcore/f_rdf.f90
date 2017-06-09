module f_rdf

use iso_fortran_env
implicit none
integer, parameter :: dp = selected_int_kind(8)
contains


subroutine dist_vec(xyz, n, box, a)
! generate pairs of distances from an xyz matrix xyz of size (n, 3)
! with periodic boundary conditions
! 29/06/16
    integer(kind=8), intent(in) :: n
    real(dp), intent(in) :: xyz(n, 3)
    real(dp), intent(out) :: a(n*(n-1)/2)
    real(dp), intent(in) :: box(3, 3)
    real(dp) :: inv_box(3, 3)
    real(dp) :: g(3), dr(3)
    real(dp) :: dist
    integer(kind=8) :: i, j, cnt

!    inv_box = box/box(1, 1)**2 ! valid only for cubic  boxes
    inv_box = 0.0
    do i = 1, 3
        inv_box(i, i) = 1.0 / box(i, i)
    enddo
 
    cnt = 1
    do i = 1, n
        do j = i+1, n
            dr(:) = xyz(i, :) - xyz(j, :)
            g = mat_vec_mul(inv_box, dr)
            g = g - nint(g)
            a(cnt) = sqrt(sum( mat_vec_mul(box, g)**2 ))
            cnt = cnt + 1
        enddo
    enddo
end subroutine


subroutine dist_vec_2mat(xyz1, n1, xyz2, n2,  box, a)
! generate pairs of mutual distances for two xyz matrices
! with periodic boundary conditions
! 29/06/16
    integer(kind=8), intent(in) :: n1, n2
    real(dp), intent(in) :: xyz1(n1, 3), xyz2(n2, 3)
    real(dp), intent(out) :: a(n1*n2)
    real(dp), intent(in) :: box(3, 3)
    real(dp) :: inv_box(3, 3)
    real(dp) :: g(3), dr(3)
    real(dp) :: dist
    integer(kind=8) :: i, j, cnt

!    inv_box = box/box(1, 1)**2 ! valid only for cubic boxes
    inv_box = 0.0
    do i = 1, size(box, 1)
        inv_box(i, i) = 1.0 / box(i, i)
    enddo
 
    cnt = 1
    do i = 1, n1
        do j = 1, n2
            dr(:) = xyz1(i, :) - xyz2(j, :)
            g = mat_vec_mul(inv_box, dr)
            g = g - nint(g)
            a(cnt) = sqrt(sum( mat_vec_mul(box, g)**2 ))
            cnt = cnt + 1
        enddo
    enddo  
end subroutine 


subroutine dist_vec_cut(xyz, n, rc, l, box, a, n_a)
! Generate pairs of distances from an xyz matrix less than cutoff rc.
! Must supply size of pair distance vector n_a externally from Python!
! Internally impossible to set if using f2py.
!
! Arguments:
! * xyz: xyz matrix of size (n, 3)
! * l: box size
! * rc: cutoff
! * box: (3, 3) matrix of cell vectors to check for nearest images
! Output:
! * a: distance vector of size (n_a)
! 05/04/16
    integer(kind=8), intent(in) :: n
    real(dp), intent(in) :: xyz(n, 3)
    real(dp), intent(in) :: rc
    real(dp), intent(in) :: l
    real(dp), intent(in) :: box(3, 3)
    integer(kind=8), intent(in) :: n_a
    real(dp), intent(out) :: a(n_a)
    real(dp) :: inv_box(3, 3)
    real(dp) :: g(3), dr(3)
    integer(kind=8) :: i, j, cnt
    real(dp) :: dist
!    real(dp), intent(out) :: a(int(2*n**2 * (rc/l)**3 * 1.1)) ! DOES NOT WORK IN f2py

!    inv_box = box/l**2   ! valid only for cubes of side l
    inv_box = 0.0
    do i = 1, size(box, 1)
        inv_box(i, i) = 1.0 / box(i, i)
    enddo
 
    cnt = 1
    do i = 1, n
        do j = i+1, n
            dr = xyz(i, :) - xyz(j, :)
            g = mat_vec_mul(inv_box, dr)
            g = g - nint(g)
            dist = sqrt(sum( mat_vec_mul(box, g)**2 ))
            if (dist < rc) then
                a(cnt) = dist
                cnt = cnt + 1
            endif
        enddo
    enddo
end subroutine


subroutine dist_vec_cut_2mat(xyz1, n1, xyz2, n2, rc, l, box, a, n_a)
! Generate pairs of distances from an xyz matrix less than cutoff rc.
!
! Arguments:
! * xyz1, xyz2: xyz matrices of size (n, 3)
! * l: box size
! * rc: cutoff
! * box: (3, 3) matrix of box vectors to check for nearest images
!
! Output:
! * a: distance vector of size (n_a)
! 05/04/16
    integer(kind=8), intent(in) :: n1, n2
    real(dp), intent(in) :: xyz1(n1, 3), xyz2(n2, 3)
    real(dp), intent(in) :: rc
    real(dp), intent(in) :: l
    real(dp), intent(in) :: box(3, 3)
    integer(kind=8), intent(in) :: n_a
    real(dp), intent(out) :: a(n_a)
    real(dp) :: inv_box(3, 3)
    real(dp) :: g(3), dr(3)
    integer(kind=8) :: i, j, cnt
    real(dp) :: dist

!    inv_box = box/l**2   ! valid only for cubes of side l
    inv_box = 0.0
    do i = 1, size(box, 1)
        inv_box(i, i) = 1.0 / box(i, i)
    enddo
 
    cnt = 1
    do i = 1, n1
        do j = i+1, n2
            dr = xyz1(i, :) - xyz2(j, :)
            g = mat_vec_mul(inv_box, dr)
            g = g - nint(g)
            dist = sqrt(sum( mat_vec_mul(box, g)**2 ))
            if (dist < rc) then
                a(cnt) = dist
                cnt = cnt + 1
            endif
        enddo
    enddo
end subroutine


subroutine histogram(a, n, n_b, hist, bins)
! Create a histogram with n_b+1 bins for vector a
! 10/01/16
    integer(kind=8), intent(in) :: n, n_b
    real(dp), intent(in) :: a(n)
    real(dp) :: dr
    integer(kind=8), intent(out) :: hist(n_b)
    real(dp), intent(out) :: bins(n_b+1)
    integer(kind=8) :: i, j

    bins(1) = minval(a)
    dr = (maxval(a) - minval(a))/n_b
    do i = 1, n_b+1
        bins(i+1) = bins(i) + dr
    enddo
    
    do j = 1, size(a)
        do i = 1, n_b
            if (a(j) > bins(i) .and. a(j) < bins(i+1)) then
                hist(i) = hist(i) + 1
                exit
            endif
        enddo
    enddo
end subroutine


subroutine pair_dist_hist(xyz, n, n_b, box, hist, bins)
! putting functions pair_dist_mat and histogram together
! 10/01/16
    integer(kind=8), intent(in) :: n, n_b
    real(dp), intent(in) :: xyz(n, 3)
    real(dp), intent(in) :: box(3, 3)
    real(dp) :: a(n*(n-1)/2)
    integer(kind=8), intent(out) :: hist(n_b)
    real(dp), intent(out) :: bins(n_b+1)

    call dist_vec(xyz, n, box, a)
    call histogram(a, n*(n-1)/2, n_b, hist, bins)
end subroutine


subroutine pair_dist_hist_2mat(xyz1, n1, xyz2, n2, n_b, box, hist, bins)
! putting functions pair_dist_mat2 and histogram together
! 10/01/16
    integer(kind=8), intent(in) :: n1, n2, n_b
    real(dp), intent(in) :: xyz1(n1, 3)
    real(dp), intent(in) :: xyz2(n2, 3)
    real(dp), intent(in) :: box(3, 3)
    real(dp) :: a(n1*n2)
    integer(kind=8), intent(out) :: hist(n_b)
    real(dp), intent(out) :: bins(n_b+1)

    call dist_vec_2mat(xyz1, n1, xyz2, n2, box, a)
    call histogram(a, n1*n2, n_b, hist, bins)
end subroutine


function mat_vec_mul(A, x) result (y)
! Simple multiplication of a (3, 3) matrix and a vector
! Should use (dp) for reals?
    real, intent(in) :: A(3, 3)
    real, intent(in) :: x(3)
    real :: y(3)
    integer :: i, j
    y = 0.0

    do i = 1, 3
        do j = 1, 3
            y(i) = y(i) + A(i, j) * x(j)
        enddo
    enddo
end function


end module


