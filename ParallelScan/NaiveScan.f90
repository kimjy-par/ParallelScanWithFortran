program prefixsum
implicit none

    integer,parameter :: num=500000000
    integer :: ARR(num), ARR_SCAN(num), ARR_GPU(num)
    integer :: i, k, ite, kte
    integer :: t1, t2, clock_rate, clock_max

    do i = 1, num
        if(MOD(i,2)==0) then 
            ARR(i) = 1
        else
            ARR(i) = 0
        end if
    end do
    call system_clock(t1, clock_rate, clock_max)
    call IntegerNaiveScan(ARR, ARR_SCAN, num)
    call system_clock(t2, clock_rate, clock_max)
    write(*,*) "Time of Naive Scan in CPU =", real(t2 - t1)/real(clock_rate)

    call system_clock(t1, clock_rate, clock_max)
    call IntegerNaiveScan_gpu(ARR, ARR_GPU, num)
    call system_clock(t2, clock_rate, clock_max)
    write(*,*) "Time of Naive Scan in GPU =", real(t2 - t1)/real(clock_rate)
    
    do i = 1, num
        if((ARR_SCAN(i)-ARR_GPU(i)) /= 0) write(*,*) "Error at ",i,"th loop!!!"
    end do
end program

subroutine IntegerNaiveScan(target_arr, result_arr, n)
implicit none
    integer, intent(in)::n
    integer, intent(in) :: target_arr(n)
    integer, intent(out) :: result_arr(n)
    real :: temp(n)
    real ::Real_ite
    integer ::i, k, ite, kte
    
    Real_ite=alog(n)/alog(2)
    ite = int(Real_ite)
    if(Real_ite-ite>0) ite = ite + 1
    kte = n
    result_arr(:) = target_arr(:)

    do i = 1, ite
    do k = 1, kte
        if(k>=2**(i-1)+1) then
            temp(k) = result_arr(k-2**(i-1)) + result_arr(k)
        else
            temp(k) = result_arr(k)
        end if
    end do
        result_arr(:) = temp(:)
    end do
end subroutine

subroutine IntegerNaiveScan_gpu(target_arr, result_arr, n)
implicit none
    integer, intent(in)::n
    integer, intent(in) :: target_arr(n)
    integer, intent(out) :: result_arr(n)
    real :: temp(n)
    real ::Real_ite
    integer ::i, k, ite, kte
    
    Real_ite=alog(n)/alog(2)
    ite = int(Real_ite)
    if(Real_ite-ite>0) ite = ite + 1
    kte = n
!$acc kernels pcopyin(target_arr) pcopyout(result_arr) &
!$acc         pcreate(temp)
!$acc loop gang vector
    do k = 1, kte
        result_arr(k) = target_arr(k)
    end do

!$acc loop seq
    do i = 1, ite
!$acc loop gang vector
    do k = 1, kte
        if(k>=2**(i-1)+1) then
            temp(k) = result_arr(k-2**(i-1)) + result_arr(k)
        else
            temp(k) = result_arr(k)
        end if
    end do
!$acc loop gang vector
    do k = 1, kte
        result_arr(k) = temp(k)
    end do
    end do
!$acc end kernels
end subroutine
