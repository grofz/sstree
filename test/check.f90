  module check_mod
    use sstree
    implicit none
    integer, parameter :: I8 = selected_int_kind(18)
    integer, parameter :: DP = kind(1.0d0)
  contains

    subroutine test3(n)
      integer, parameter :: cycles=120
      integer, intent(in) :: n
      type(point_t) :: p0, p1
      type(point_t), allocatable :: p(:)
      integer, allocatable :: ind(:)
      integer :: i, j
      type(sstree_t) :: tree
      logical :: isvalid
      real, dimension(CYCLES) :: s0, s1, e0, e1

      p0 % x = 0.0
      p1 % x = 10.0
      call tree % insert(p0)
      goto 100
      p = random_points(n/10, p0, p1)
      ind = random_order(n/10)
      do i= 1, n/10 
        call tree % insert(p(ind(i)))
      enddo
      isvalid = tree % isvalid()
      if (.not. isvalid) error stop 'invalidate during insert'
print *, ' initial 10% of nodes'

      100 continue
      do j = 1, CYCLES
        p = random_points(n, p0, p1)
        ind = random_order(n)

        call cpu_time(s0(j))
        do i= 1, n 
!  print '(a,3(f5.2,1x))', 'insert =',p(ind(i))%x
          call tree % insert(p(ind(i)))
          !isvalid = tree % isvalid()
          !if (.not. isvalid) error stop 'invalidate during insert'
!  print *
        enddo
        call cpu_time(e0(j))
          isvalid = tree % isvalid()
          if (.not. isvalid) error stop 'invalidate during insert'

        print *, 'delete'
        call cpu_time(s1(j))
        do i = 1, n  
!  print '(a,3(f5.2,1x))', 'delete =',p(i)%x
          call tree % delete(p(i))
          !isvalid = tree % isvalid()
          !if (.not. isvalid) error stop 'invalidate during delete'
          !print *

          !if (mod(i,n/20) /= 0) cycle
          isvalid = tree % isvalid()
          if (.not. isvalid) error stop 'invalidate during insert'
        enddo
        call cpu_time(e1(j))
        print *, 'end of cycle ',j
      enddo

      do j=1,CYCLES
        print *, e0(j)-s0(j), e1(j)-s1(j)
      enddo




    end subroutine test3

    subroutine test1()
      type(sstree_t) :: tree
      type(point_t) :: p(14)
      logical :: isvalid
      integer :: i

      p(1) = point_t([1.0, 0.0, 0.0])
      p(2) = point_t([1.0, 1.0, 0.0])
      p(3) = point_t([0.0, 0.0, 0.0])
      p(4) = point_t([0.5, 0.5, 0.0])
      p(5) = point_t([0.1, 0.1, 0.0])
      p(6) = point_t([0.9, 0.1, 0.0])
      p(7) = point_t([0.9, 0.9, 0.0])
      p(8) = point_t([1.0, 0.1, 0.0])
      p(9) = point_t([1.0, 0.9, 0.0])
      p(10) = point_t([0.2, 0.2, 0.0])
      p(11) = point_t([0.2, 0.3, 0.0])
      p(12) = point_t([0.2, 0.4, 0.0])
      p(13) = point_t([0.9, 0.7, 0.0])
      p(14) = point_t([0.1, 0.7, 0.0])

      isvalid = tree % isvalid()
      do i=1, 14 
   print '(a,3(f5.2,1x))', 'insert =',p(i)%x
        call tree % insert(p(i))
        isvalid = tree % isvalid()
   print *
      enddo
      !stop

      print *, 'delete'
      do i=1, 14  
        call tree % delete(p(i))
        isvalid = tree % isvalid()
        print *
      enddo

    end subroutine



    subroutine test2(n, inbetween_validation, time)
      integer, intent(in) :: n
      type(sstree_t) :: tree
      type(point_t) :: p
      logical :: isvalid, inbetween_validation
      integer, allocatable :: seed(:)
      integer :: i
      real(DP) :: time, tstart, tend

      goto 100
      call random_seed(size=i)
      print *, 'seed size =', i
      allocate(seed(i))
      call random_seed(get=seed)
      print *, 'old seed =',seed
      seed = 12345
      call random_seed(put=seed)

      100 continue
      isvalid = tree % isvalid()
      call cpu_time(tstart)
      do i=1, n 
        call random_number(p % x)
  !print '(a,3(f5.2,1x))', 'insert =',p%x
        call tree % insert(p)
        if (inbetween_validation) then
          isvalid = tree % isvalid()
          if (.not. isvalid) error stop 'tree invalidated'
        endif
  !print *
      enddo
      call cpu_time(tend)
      isvalid = tree % isvalid()
      if (.not. isvalid) error stop 'resulting tree not valid'
      time = tend - tstart
      print '("Insert took ",en12.3," s ",en12.3," s/item")', &
          time, time / n

    end subroutine test2



    function random_points(n,p0,p1) result(points)
      integer, intent(in) :: n
      type(point_t), intent(in) :: p0, p1
      type(point_t), allocatable :: points(:)

      integer :: i
      allocate(points(n))
      do i=1, size(points(1) % x)
        call random_number(points % x(i))
        points % x(i) = p1%x(i) * points % x(i) &
                      + p0%x(i) * (1.0 - points % x(i))
      enddo
    end function random_points



    function random_order(n) result(arr)
      integer, intent(in) :: n
      integer, allocatable :: arr(:)

      integer :: i, isel, tmp
      real    :: rnd

      allocate(arr(n))
      do i = 1, n
        arr(i) = i
      enddo

      do i = 1, n-1
        call random_number(rnd)
        isel = int((n-i+1)*rnd) ! isel = (0, n-i)
        tmp = arr(i)
        arr(i) = arr(i+isel)
        arr(i+isel)=tmp
      enddo
    end function random_order

  end module check_mod

  program check
    use check_mod
    implicit none
    integer :: i
    real(DP) :: time
    !call test1()
    call test3(900)
    stop

    do i=1,1
      call test2(20000, .true., time)
    enddo
    !call test2(1000000, .false., time)
    print *, 'Tests finished!'
  end program check
