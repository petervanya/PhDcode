module f_rdf2
! Create histograms for RDFs of one or two particle types
! 08/12/17

implicit none
contains

subroutine rdf_hist(xyz, dr, box, nb, h)
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: xyz(:, :)
    real(dp), intent(in) :: dr
    integer(8), intent(in) :: nb
    real(dp), intent(in) :: box(3, 3)
    integer(8), intent(out) :: h(nb)
    real(dp) :: inv_box(3, 3), box_diag(3), g(3), r(3), d
    integer(8) :: i, j

    h(:) = 0
    inv_box(:, :) = 0.0
    do i = 1, 3
        inv_box(i, i) = 1.0 / box(i, i)
        box_diag(i) = box(i, i)
    end do
 
    do i = 1, size(xyz, 1)
        do j = i+1, size(xyz, 1)
            r(:) = xyz(i, :) - xyz(j, :)
            g = matmul(inv_box, r)
            g = g - nint(g)
            d = sqrt(sum(matmul(box, g)**2))
            h(int(d / dr) + 1) = h(int(d / dr) + 1) + 1
        enddo
    enddo
end subroutine

subroutine rdf_hist_2_types(xyz1, xyz2, dr, box, nb, h)
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: xyz1(:, :), xyz2(:, :)
    real(dp), intent(in) :: dr
    integer(8), intent(in) :: nb
    real(dp), intent(in) :: box(3, 3)
    integer(8), intent(out) :: h(nb)
    real(dp) :: inv_box(3, 3), box_diag(3), g(3), r(3), d
    integer(8) :: i, j

    h(:) = 0
    inv_box(:, :) = 0.0
    do i = 1, 3
        inv_box(i, i) = 1.0 / box(i, i)
        box_diag(i) = box(i, i)
    end do
 
    do i = 1, size(xyz1, 1)
        do j = 1, size(xyz2, 1)
            r(:) = xyz1(i, :) - xyz2(j, :)
            g = matmul(inv_box, r)
            g = g - nint(g)
            d = sqrt(sum(matmul(box, g)**2))
            h(int(d / dr) + 1) = h(int(d / dr) + 1) + 1
        enddo
    enddo  
end subroutine 

end module
