module mat_ops
implicit none
contains

subroutine get_pair_dist(mat, n, arr)
! generate pairs of distances from an xyz matrix mat of size (n, 3)
! 22/09/15
    integer, parameter :: k = selected_int_kind(16)
    integer(kind=k), intent(in) :: n
    real(8), intent(in), dimension(n, 3) :: mat
    real(8), intent(out), dimension(n*(n-1)/2) :: arr
    integer :: i, j, cnt
    ! f2py depend(n) mat, arr

    cnt = 1
    do i = 1, n
        do j = i+1, n
            arr(cnt) = sqrt(sum((mat(i, :) - mat(j, :))**2))
            cnt = cnt + 1
        enddo
    enddo
end subroutine


subroutine get_pair_dist2(mat, N, arr)
! generate pairs of distances together with the particle types
! 11/10/15
    integer, intent(in) :: N
    real, intent(in), dimension(N, 4) :: mat
    real, intent(out), dimension(N*(N-1)/2, 3) :: arr
    integer :: i, j, cnt

    cnt = 1
    do i = 1, N
        do j = i+1, N
            arr(cnt, 1) = mat(i, 1)
            arr(cnt, 2) = mat(j, 1)
            arr(cnt, 3) = sqrt((mat(i, 2)-mat(j, 2))**2 + &
                          (mat(i, 3)-mat(j, 3))**2 + &
                          (mat(i, 4)-mat(j, 4))**2)
            cnt = cnt + 1
        enddo
    enddo
end subroutine

subroutine get_local_op(xyz, N, rc, phi)
! Calculate and return order parameter (according to Goyal PhD thesis)
! 12/10/15
    integer, intent(in) :: N
    real, intent(in) :: rc
    real, intent(in), dimension(N, 4) :: xyz
    real, intent(out) :: phi
    integer, dimension(N) :: n1, n2
    real :: dist
    integer :: i, j
    
    n1 = 0
    n2 = 0
    phi = 0.0

    do i = 1, N
        do j = i+1, N
            dist = sqrt((xyz(i, 2)-xyz(j, 2))**2 + &
                        (xyz(i, 3)-xyz(j, 3))**2 + &
                        (xyz(i, 4)-xyz(j, 4))**2)
            if (dist < rc) then
                if (xyz(j, 1) == 1) then   ! neighbour of i-th particle is type 1
                    n1(i) = n1(i) + 1
                endif
                if (xyz(j, 1) == 2) then   ! neighbour of i-th particle is type 2
                    n2(i) = n2(i) + 1
                endif
                if (xyz(i, 1) == 1) then   ! vice versa
                    n1(j) = n1(j) + 1
                endif
                if (xyz(i, 1) == 2) then   ! vice versa
                    n2(j) = n2(j) + 1
                endif
            endif
        enddo
    enddo
    
    do i = 1, N
        phi = phi + float(n1(i) - n2(i))**2/(n1(i) + n2(i))**2
!        if (n1(i) == 0) then 
!            print *, "n1", i
!        endif
!        if (n2(i) == 0) then
!            print *, "n2", i
!        endif 
    enddo
!    print *, phi
    phi = phi/N
end subroutine


subroutine get_local_op2(xyz, N, rc, n1, n2)
! Return two nearest neighbour occupation arrays n1, n2 to later calculate 
! the order parameter (according to Goyal PhD thesis)
! 12/10/15
    integer, intent(in) :: N
    real, intent(in) :: rc
    real, intent(in), dimension(N, 4) :: xyz
    integer, intent(out), dimension(N) :: n1, n2
    real :: dist, phi
    integer :: i, j
    
    n1 = 0
    n2 = 0
    phi = 0.0

    do i = 1, N
        do j = i+1, N
            dist = sqrt((xyz(i, 2)-xyz(j, 2))**2 + &
                        (xyz(i, 3)-xyz(j, 3))**2 + &
                        (xyz(i, 4)-xyz(j, 4))**2)
            if (dist < rc) then
                if (xyz(j, 1) == 1) then   ! neighbour of i-th particle is type 1
                    n1(i) = n1(i) + 1
                endif
                if (xyz(j, 1) == 2) then   ! neighbour of i-th particle is type 2
                    n2(i) = n2(i) + 1
                endif
                if (xyz(i, 1) == 1) then   ! vice versa
                    n1(j) = n1(j) + 1
                endif
                if (xyz(i, 1) == 2) then   ! vice versa
                    n2(j) = n2(j) + 1
                endif
            endif
        enddo
    enddo
end subroutine


function get_n_bins(a, n, dr) result (n_bins)
    integer, intent(in) :: n
    real(8), intent(in) :: a(n)
    real(8), intent(in) :: dr
    integer :: n_bins
    n_bins = int((maxval(a) - minval(a))/dr) + 1
end function


pure function n_bins(a, n, dr)
! Return the number of bins based on min and max value of a and bin size dr
    integer, intent(in) :: n
    real(8), intent(in) :: a(n)
    real(8), intent(in) :: dr
    integer :: n_bins
    n_bins = int((maxval(a) - minval(a))/dr) + 1
end function


subroutine create_hist(a, n, n_b, hist, bins)
! Create histogram with bins of size dr, version for f2py
    integer, intent(in) :: n
    integer, intent(in) :: n_b
    real(8), intent(in) :: a(n)
    real(8) :: dr
    real(8), intent(out) :: bins(n_b+1)
! bins(int(sqrt(real(n)))+2) !(n_bins(a, n, dr)+1)
    integer, intent(out) :: hist(n_b)
! hist(int(sqrt(real(n)))+1) !(n_bins(a, n, dr))
    integer :: i, j
    hist = 0
    
    bins(1) = minval(a)
    dr = (maxval(a) - minval(a))/n_b
    do i = 1, n_b
        bins(i+1) = bins(i) + dr
    enddo
    
    ! Create the histogram
    do j = 1, size(a)
        do i = 1, n_b
            if (a(j) > bins(i) .and. a(j) < bins(i+1)) then ! CHECK
                hist(i) = hist(i) + 1
                exit
            endif
        enddo
    enddo
end subroutine


!subroutine create_hist_gnu(a, n, dr, bins, hist)
!! Create histogram with bins of size dr, version for gfortran
!    integer, intent(in) :: n
!    real(8), intent(in) :: a(n)
!    real(8), intent(in) :: dr
!    integer, intent(out), allocatable :: hist(:)
!    real(8), intent(out), allocatable :: bins(:)
!    integer :: i, j
!    
!    n_b = n_bins(a, n, dr)
!    allocate(bins(n_b+1))
!    allocate(hist(n_b))
!    hist = 0
!    
!    ! Create bins
!    bins(1) = minval(a)
!    do i = 1, n_b
!        bins(i+1) = bins(i) + dr
!    enddo
!    print *, bins
!    print *, hist
!    
!    ! Create the histogram
!    do j = 1, n
!        do i = 1, n_b
!            if (a(j) > bins(i) .and. a(j) < bins(i+1)) then ! CHECK
!                hist(i) = hist(i) + 1
!                exit
!            endif
!        enddo
!    enddo
!end subroutine
end module


