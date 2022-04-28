  module check_mod
    use sstree
    use point_mod, only : point_t, sphere_t, rectangle_t
    implicit none
  contains

    subroutine test1()
      type(sstree_t) :: tree
      type(point_t) :: p(14), tar, fou
      type(rectangle_t) :: rect
      type(sphere_t) :: sph
      type(point_t), allocatable :: psel(:)
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

      !tar = point_t([0.92, 0.92, 0.0])
      !fou = tree % nearestNeighbor(tar)
      !print '(3(f5.2,1x))', fou % x

      rect % lcor % x = [0.0, 0.0, 0.0]
      rect % rcor % x = [0.5, 0.5, 0.0]
      sph % center % x = [0.0, 0.0, 0.0]
      sph % radius = 0.14143
      print *,'pointsWithinRegion...'
      !psel = tree % pointsWithinRegion(rect)
      psel = tree % pointsWithinRegion(sph)
      print *,'...ok'

      do i = 1, size(psel)
        print '(3(f5.2,1x))', psel(i)%x
      enddo
      if (size(psel)==0) print *, 'no points in the region' 

      !call delete
    contains
      subroutine delete
        print *, 'delete'
        do i=1, 14 
          call tree % delete(p(i))
          isvalid = tree % isvalid()
          print *
        enddo
      end subroutine
    end subroutine test1

  end module check_mod



  program check
    use check_mod, only : test1
    use test_sstree_mod
    use utest_mod, only : utest_t
    implicit none
    type(utest_t) :: utest

    utest = utest_t()
    call test1()
    call testInsertDelete(utest, 10000, 10, 3)
    call testNNSearch(utest, 10000, 100, 10, 6)

    call utest % summarize()
    print *, 'Tests finished!'
  end program check
