  module test_sstree_mod
    use sstree
    use utest_mod
    implicit none
    private

    public testInsertDelete, testNNSearch
 public random_points, random_order, bf_nearestNeighbor
    integer, parameter :: DP = kind(1.0d0)

  contains

    subroutine testInsertDelete(utest, n, ntimes, brfac)
      type(utest_t), intent(inout) :: utest
      integer, intent(in) :: n, ntimes
      integer, intent(in), optional ::  brfac
!
! Test "insert" and "delete" operations.
! Load the tree with "n" random points, then delete them.
! Repeat "ntimes" times. 
!
      type(sstree_t) :: tree
      type(point_t) :: p0, p1
      type(point_t), allocatable :: p(:)
      integer, allocatable :: ind(:)
      integer :: i, j, nmore
      logical :: isvalid
      real :: tstart, tend
      real, dimension(ntimes) :: t_insert, t_delete

      ! New tree
      if (present(brfac)) then
        tree = sstree_t(brfac)
      else
        tree = sstree_t(2) ! default m0=2, m1=3
      endif

      ! Points generated within the box [P0 | P1]
      p0 % x = 0.0
      p1 % x = 10.0

      ! Any additional points that remains in tree?
      nmore = max(1, n/100)

      p = random_points(nmore, p0, p1)
      do i= 1, nmore 
        call tree % insert(p(i))
      enddo
      isvalid = tree % isvalid()
      call utest % assert_eq(isvalid, .true., 'pre-loading ok?')
      if (.not. isvalid) return

      ! Main loop
      do j = 1, ntimes
        p = random_points(n, p0, p1)
        ind = random_order(n)

        print '("Loading tree (n=",i0,") ...")', n
        call cpu_time(tstart)
        do i= 1, n 
          call tree % insert(p(ind(i)))
        enddo
        call cpu_time(tend)
        t_insert(j) = tend - tstart
        isvalid = tree % isvalid()
        if (j==1) &
          call utest % assert_eq(isvalid, .true., 'insert ok?')
        if (.not. isvalid) exit

        print '(a)', 'Unloading tree...'
        t_delete(j) = 0.0
        call cpu_time(tstart)
        do i = 1, n  
          call tree % delete(p(i))

          ! validate 10x times and for last 40 nodes
          if (mod(i,n/10) == 0 .or. n-i < 40) then
            call cpu_time(tend)
            t_delete(j) = t_delete(j) + (tend - tstart)
            isvalid = tree % isvalid()
            call cpu_time(tstart)
            if (.not. isvalid) exit
          endif
        enddo
        call cpu_time(tend)
        t_delete(j) = t_delete(j) + (tend - tstart)
        if (j==1) &
          call utest % assert_eq(isvalid, .true., 'delete ok?')

        print *, 'end of cycle ',j
        if (.not. isvalid) exit
      enddo

      do j=1, ntimes
        print *, t_insert(j), t_delete(j)
      enddo

      call utest % assert_eq(isvalid, .true., 'insert-delete cycles ok?')
    end subroutine testInsertDelete



    subroutine testNNSearch(utest, n, m, ntimes, brfac)
      type(utest_t), intent(inout) :: utest
      integer, intent(in) :: n, m, ntimes
      integer, intent(in), optional :: brfac
!
! Test "nearestNeighbor" search. Generate and load "n" points to
! tree, generate testing sample of "m" points. Compare search
! using tree and brute-force search. Repeat "ntimes". In every
! cycle, replace "n/REPL_FACTOR" points in the tree.
!
      integer, parameter :: REPL_FACTOR = 5
      type(sstree_t) :: tree
      type(point_t), allocatable :: p(:), q(:), f(:), ff(:), prepl(:)
      type(point_t) :: p0, p1
      integer, allocatable :: ind(:)
      logical :: isvalid
      integer :: i, itime, nrepl
      real(DP) :: tstart, tend
      real(DP) :: t(ntimes), tt(ntimes)

      ! New tree
      if (present(brfac)) then
        tree = sstree_t(brfac)
      else
        tree = sstree_t(2)
      endif

      ! "n" points generated within the box [P0 | P1] and loaded to tree
      p0 % x = 0.0
      p1 % x = 10.0
      p = random_points(n, p0, p1)
      ind = random_order(n)
      print *, 'loading tree...'
      do i = 1, n
        call tree % insert(p(ind(i)))
      enddo
      isvalid = tree % isvalid()
      if (.not. isvalid) goto 100

      t = 0.0  ! time elapsed using sstree search
      tt = 0.0 ! time elapsed using brute-force search

      ! Main loop
      do itime = 1, ntimes

        ! Generate sample of "m" points to be searched for
        q = random_points(m, p0, p1)
        allocate(f(m), ff(m))

        ! SS tree search
        call cpu_time(tstart)
        do i = 1, m
          f(i) = tree % nearestNeighbor(q(i))
        enddo
        call cpu_time(tend)
        t(itime) = tend - tstart

        ! Brute force search
        call cpu_time(tstart)
        do i = 1, m
          ff(i) = bf_nearestNeighbor(p,q(i))
        enddo
        call cpu_time(tend)
        tt(itime) = tend - tstart

        ! Both methods found same points?
        isvalid = all(f == ff)
        if (itime < 3) &
            call utest % assert_eq(isvalid, .true., 'NN search - same points found?')
        if (.not. isvalid) exit
        deallocate(f, ff)

        ! Delete some points from tree and replace them by a new set
        if (itime == ntimes) exit
        nrepl = n / REPL_FACTOR
        prepl = random_points(nrepl, p0, p1)
        ind = random_order(nrepl)
        do i = 1, nrepl
          call tree % delete(p(i))
          call tree % insert(prepl(ind(i)))
        enddo
        isvalid = tree % isvalid()
        if (.not. isvalid) goto 100
        p(1:nrepl) = prepl
        print *, 'end of cycle =', itime
      enddo

      print *
      print '(a)','Nearest neighbor searching time'
      print '(a)','     SS tree  brute force     speed-up'
      print '(a)','--------------------------------------'
      do i = 1, ntimes
        print 301, t(i), tt(i), tt(i)/t(i)
        !301 format (2(en12.3,1x),f12.3,1x)
        301 format (2(f12.3,1x),f12.3,1x)
      enddo
      print *

      100 call utest % assert_eq(isvalid, .true., 'NN search cycles ok?')
    end subroutine testNNSearch



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



    function bf_nearestNeighbor(p, q) result(nn)
      type(point_t) :: nn
      type(point_t), intent(in) :: p(:), q

      integer  :: i
      real(DP) :: nn_dist, dist

      nn_dist = huge(nn_dist)
      do i = 1, size(p)
        dist = p(i) % distance(q)
        if (dist >= nn_dist) cycle
        nn = p(i)
        nn_dist = dist
      enddo
    end function bf_nearestNeighbor

  end module test_sstree_mod
