  module check_mod
    use sstree
    implicit none
    integer, parameter :: I8 = selected_int_kind(18)
    integer, parameter :: DP = kind(1.0d0)
  contains

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

  end module check_mod

  program check
    use check_mod
    implicit none
    integer :: i
    real(DP) :: time
    !call test1()

    do i=1,1
      call test2(200, .true., time)
    enddo
    !call test2(1000000, .false., time)
    print *, 'Tests finished!'
  end program check
