  module sstree
    use iso_fortran_env, only : output_unit
    implicit none
    private
    public sstree_t, point_t

    integer, parameter :: DIM=3, M0=1, M1=3
    integer, parameter :: WP = kind(1.0d0)
    !integer, parameter :: MSG_MEMORY = output_unit
    integer, parameter :: MSG_MEMORY = 0

    type point_t
      real(WP) :: x(DIM)
    contains
      procedure :: distance => point_distance
      procedure :: isin => point_isin
      procedure, private :: point_equals
      generic :: operator(==) => point_equals
    end type

    type ssnode_ptr
      type(ssnode_t), pointer :: p => null()
    end type

    type ssnode_t
      type(point_t)    :: centroid
      real(WP)         :: radius
      type(ssnode_ptr) :: children(M1+1)
      type(point_t)    :: points(M1+1)
      logical          :: isleaf ! ... or node with children
      integer          :: n = 0
    contains
      procedure :: intersectsPoint => ssnode_intersectsPoint
      procedure :: findClosestChild => ssnode_findClosestChild
      procedure :: updateBoundingEnvelope => ssnode_updateBoundingEnvelope
      procedure :: directionOfMaxVariance => ssnode_directionOfMaxVariance
      procedure :: findSplitIndex => ssnode_findSplitIndex
      procedure :: getCentroids => ssnode_getCentroids
      procedure :: minVarianceSplit => ssnode_minVarianceSplit
      procedure :: getRadii => ssnode_getRadii
      procedure :: addPoint => ssnode_addPoint
      procedure :: addChild => ssnode_addChild
      procedure :: deleteChild => ssnode_deleteChild
      procedure :: split => ssnode_split
    end type

    interface ssnode_t
      module procedure new_node
      module procedure new_leaf
    end interface

    type sstree_t
      type(ssnode_t), pointer :: root => null()
      integer :: m0, m1, dim
    contains
      procedure :: insert => sstree_insert
      procedure :: isvalid => sstree_isvalid
      final :: sstree_finalize
    end type

    interface sstree_t
      module procedure new_tree
    end interface

  contains

    function new_node(children)
      type(ssnode_t), pointer :: new_node
      type(ssnode_ptr), intent(in) :: children(:)

      allocate(new_node)
      new_node % n = size(children)
      new_node % isleaf = .false.
      new_node % children(1: new_node % n) = children
      call new_node % updateBoundingEnvelope()
      if (MSG_MEMORY /= 0) then
        write(MSG_MEMORY, '(a,i0,a,3(en12.3,1x),a,en12.3)') &
        '(newnode) children=', new_node % n, &
        '  XC=',new_node%centroid%x,'  R=',new_node%radius
      endif
    end function new_node



    function new_leaf(points)
      type(ssnode_t), pointer :: new_leaf
      type(point_t), intent(in) :: points(:)

      allocate(new_leaf)
      new_leaf % n = size(points)
      new_leaf % isleaf = .true.
      new_leaf % points(1: new_leaf % n) = points
      call new_leaf % updateBoundingEnvelope()
      if (MSG_MEMORY /= 0) then
        write(MSG_MEMORY, '(a,i0,a,3(en12.3,1x),a,en12.3)') &
        '(newleaf) children=', new_leaf % n, &
        '  XC=',new_leaf%centroid%x,'  R=',new_leaf%radius
      endif
    end function new_leaf



    function new_tree(adim, am0, am1)
      type(sstree_t) :: new_tree
      integer, intent(in) :: adim, am0, am1

      new_tree % m0 = am0
      new_tree % m1 = am1
      new_tree % dim = adim
! TODO not finished
    end function



! Listing 10.3
    recursive function search(node, target) result(found_node)
      type(ssnode_t), pointer :: found_node
      type(ssnode_t), intent(in), pointer :: node
      type(point_t), intent(in) :: target
!
! Return leaf that contains the target point, return null otherwise.
!
      integer :: i

      if (.not. associated(node)) &
          error stop 'search - unassociated node'

      if (node % isleaf) then
        do i = 1, node % n
          if (node % points(i) == target) then
            found_node => node
            return
          endif
        enddo     
      else
        do i = 1, node % n
          if (.not. node % children(i)%p % intersectsPoint(target)) cycle 
          found_node => search(node % children(i)%p, target)
          if (associated(found_node)) return
        enddo
      endif
      found_node => null()
    end function search



! Listing 10.4
    pure logical function ssnode_intersectsPoint(this, point) result(res)
      class(ssnode_t), intent(in) :: this
      class(point_t), intent(in) :: point

      res = this % centroid % distance(point) <= this % radius
    end function ssnode_intersectsPoint



    elemental function point_distance(p0, p1) result(dis)
      real(WP) :: dis
      class(point_t), intent(in) :: p0, p1
!
! Euclidian distance between points P0 and P1
!
      dis = sum((p0 % x - p1 % x)**2)
      dis = sqrt(dis)
    end function point_distance



    pure logical function point_equals(p0, p1)
      class(point_t), intent(in) :: p0, p1
      point_equals = all(p0 % x == p1 % x)
    end function point_equals



! Listing 10.5
    recursive function searchParentLeaf(node, target) result(found_node)
      type(ssnode_t), pointer :: found_node
      type(ssnode_t), pointer, intent(in) :: node
      type(point_t), intent(in) :: target
!
! Return closest leaf to a target point
!
      type(ssnode_t), pointer :: child

      if (.not. associated(node)) &
        error stop 'searchParentLeaf - unassociated node'

      if (node % isleaf) then
        found_node => node
      else
        child => node % findClosestChild(target)
        found_node => searchParentLeaf(child, target)
      endif
    end function searchParentLeaf



 ! Listing 10.8
    function ssnode_findClosestChild(this, point) result(res)
      class(ssnode_t), pointer :: res
      class(ssnode_t), intent(in) :: this
      type(point_t), intent(in) :: point
 !
 ! Return a child of a current node whose distance to point is minimal
 !
      real(WP) :: mindis, dis
      integer :: i

      if (this % isleaf) &
          error stop 'findClosestChild - node is leaf'
      mindis = huge(mindis)
      res => null()
      do i = 1, this % n
        dis = this % children(i) % p % centroid % distance(point)
        if (dis >= mindis) cycle
        mindis = dis
        res => this % children(i) % p
      enddo

      if (.not. associated(res)) & ! TODO can be removed after debug
          error stop 'findClosestChild found nothing'
    end function ssnode_findClosestChild



 ! Listing 10.6
    recursive function insert(node, point) result(res)
      type(ssnode_ptr) :: res(2)
      type(ssnode_t), pointer :: node
      type(point_t), intent(in) :: point
!
! Insert point to node's subtree. Returns null (if all ok) or
! a pair of nodes resulting from the split.
!
      class(ssnode_t), pointer :: child
      type(ssnode_ptr) :: res0(2)

      if (.not. associated(node)) &
          error stop 'insert - unassociated node'

      !res(1) % p => null() ! TODO should be null implicitly
      !res(2) % p => null()

      if (node % isleaf) then
        if (point % isin(node % points(1:node % n))) then
 print *, 'A =', point % x
 print *, 'Bx=', node % points(1:node %n) % x(1)
 print *, 'By=', node % points(1:node %n) % x(2)
          print *, 'insert - point is already in tree'
          return
        endif
        call node % addPoint(point)
        call node % updateBoundingEnvelope()
        if (node % n <= M1) return
      else
        child => node % findClosestChild(point)
        res0 = insert(child, point)
        if (.not. associated(res0(1) % p)) then
          call node % updateBoundingEnvelope()
          return
        else
          call node % deleteChild(child)
          call node % addChild(res0(1) % p)
          call node % addChild(res0(2) % p)
          call node % updateBoundingEnvelope()
          if (node % n <= M1) return
        endif
      endif
      res = node % split()
    end function insert



 ! Listing 10.7
    subroutine sstree_insert(this, point)
      class(sstree_t) :: this
      type(point_t), intent(in) :: point
!
! Insert point to SS tree. Insert first point straight away, use recursive
! "ïnsert" function otherwise.
! 
      type(ssnode_ptr) :: new(2)

      if (.not. associated(this % root)) then
        this % root => ssnode_t([point])
        call this % root % updateBoundingEnvelope() !TODO move to constructor?
      else
        new = insert(this % root, point)
        if (associated(new(1) % p)) then
          if (MSG_MEMORY /= 0) then
            write(MSG_MEMORY, '(a,i0,a,l1)') &
            '(insert dealloc root) nold=', this%root%n, &
            ' isleaf=', this%root%isleaf
          endif
          deallocate(this % root)
          this % root => ssnode_t([new(1), new(2)])
        endif
      endif
    end subroutine sstree_insert



    pure logical function point_isin(this, arr)
      class(point_t), intent(in) :: this
      class(point_t), intent(in) :: arr(:)

      integer :: i
      
      point_isin = .false.
      do i = 1, size(arr)
        if (arr(i) == this) then
          point_isin = .true.
          exit
        endif
      enddo
    end function



    subroutine ssnode_addPoint(this, point)
      class(ssnode_t), intent(inout) :: this
      class(point_t), intent(in) :: point

      if (.not. this % isleaf) &
          error stop 'addPoint - is not leaf'
      if (this % n > M1) &
          error stop 'addPoint - size oveflows'
      if (point % isin(this % points(1:this % n))) &
          error stop 'addPoint - duplicate point'

      this % n = this % n + 1
      this % points(this % n) = point
    end subroutine



! Listing 10.10.
    subroutine ssnode_updateBoundingEnvelope(this)
      class(ssnode_t), intent(inout) :: this

      type(point_t), allocatable :: points(:)
      real(WP), allocatable :: radii(:)
      real(WP) :: dis
      integer :: i

      points = this % getCentroids()
      do i = 1, DIM
        this % centroid % x(i) = mean(points, i)
      enddo

      this % radius = 0.0
      if (this % isleaf) then
        do i = 1, this % n
          dis = points(i) % distance(this % centroid)
          if (dis > this % radius) this % radius = dis
        enddo
      else
        radii = this % getRadii()
        do i = 1, this % n
          dis = points(i) % distance(this % centroid) + radii(i)
          if (dis > this % radius) this % radius = dis
        enddo
      endif
    end subroutine ssnode_updateBoundingEnvelope



    subroutine ssnode_addChild(this, child)
      class(ssnode_t), intent(inout) :: this
      class(ssnode_t), pointer, intent(in) :: child

      integer :: i

      if (this % isleaf) &
          error stop 'addChild - is not internal node'
      if (this % n > M1) &
          error stop 'addChild - size oveflows'
      do i = 1, this % n
        if (associated(this % children(i) % p, child)) &
            error stop 'addChild - duplicate node'
      enddo

      this % n = this % n + 1
      this % children(this % n) % p => child
    end subroutine



    subroutine ssnode_deleteChild(this, child)
      class(ssnode_t), intent(inout) :: this
      class(ssnode_t), pointer, intent(inout) :: child

      integer :: i, i0

      if (this % isleaf) &
          error stop 'deleteChild - is not internal node'

      i0 = 0
      do i = 1, this % n
        if (associated(this % children(i) % p, child)) then
          i0 = i
        else
          if (i0 /= 0) this % children(i-1) % p => this % children(i) % p
        endif
      enddo
  
      if (i0 == 0) &
          error stop 'deleteChild - to be deleted node not found'

      this % children(this % n) % p => null()
      this % n = this % n - 1

      if (MSG_MEMORY /= 0) then
        write(MSG_MEMORY, '(a,i0,a,l1)') &
        '(deletechild) nold=', child%n, ' isleaf=', child%isleaf
      endif
      deallocate(child)
    end subroutine



!Listing 10.11
    function ssnode_split(this) result(split)
      type(ssnode_ptr) :: split(2)
      class(ssnode_t), intent(inout) :: this

      integer :: ind

      call this % findSplitIndex(ind)
      if (this % isleaf) then
        split(1) % p => ssnode_t(this % points(1:ind-1))
        split(2) % p => ssnode_t(this % points(ind:this % n))
      else
        split(1) % p => ssnode_t(this % children(1:ind-1))
        split(2) % p => ssnode_t(this % children(ind: this % n))
      endif
    end function



!Listing 10.12
    subroutine ssnode_findSplitIndex(this, split_index)
      class(ssnode_t), intent(inout) :: this
      integer, intent(out) :: split_index
!
! Prepare for splitting the node. Reorder points/children and return
! index where to split.
!
      integer :: dir, i, j
      real(WP) :: key
      type(point_t) :: tmp_point
      type(ssnode_t), pointer :: tmp_node

      dir = this % directionOfMaxVariance()

      ! insertion sort for points or node pointers
      if (this % isleaf) then
        do i = 1, this % n
          tmp_point = this % points(i)
          key = tmp_point % x(dir)
          j = i - 1
          do while (j >= 1)
            if ( this % points(j) % x(dir) <= key) exit
            this % points(j + 1) = this % points(j)
            j = j - 1
          enddo
          this % points(j + 1) = tmp_point
        enddo
      else
        do i = 1, this % n
          tmp_node => this % children(i)%p
          key = tmp_node % centroid % x(dir)
          j = i - 1
          do while (j >= 1)
            if (this%children(j)%p%centroid%x(dir) <= key) exit
            this % children(j + 1)%p => this % children(j)%p
            j = j - 1
          enddo
          this % children(j + 1)%p => tmp_node
        enddo
      endif

      split_index = this % minVarianceSplit(dir)
    end subroutine ssnode_findSplitIndex



!Listing 10.14
    pure function ssnode_minVarianceSplit(this, dir) result(split_index)
      integer :: split_index
      class(ssnode_t), intent(in) :: this
      integer, intent(in) :: dir
! 
! Get spliting index - index where the second half begins.
! Assuming that points/children arrays are sorted along dir's dimension
!
      type(point_t), allocatable :: points(:)
      real(WP) :: minvar, var
      integer :: i
      
      points = this % getCentroids()
      minvar = huge(minvar)
      split_index = M0 + 1
      do i = M0 + 1, this % n - M0 + 1
        var  = variance(points(1:i-1), dir) &
            &+ variance(points(i:this % n), dir)
        if (var < minvar) then
          minvar = var
          split_index = i
        endif
      enddo
    end function ssnode_minVarianceSplit



! Listing 10.9
    integer function ssnode_directionOfMaxVariance(this) result(dir)
      class(ssnode_t), intent(in) :: this

      type(point_t), allocatable :: points(:) 
      real(WP) :: maxvar, var
      integer :: i

      dir = 1
      maxvar = 0.0
      points = this % getCentroids()
      do i = 1, DIM
        var = variance(points, i)
        if (var > maxvar) then
          maxvar = var
          dir = i
        endif
      enddo
    end function ssnode_directionOfMaxVariance



! Listing 10.13
    pure function ssnode_getcentroids(this) result(centroids)
      class(ssnode_t), intent(in) :: this
      type(point_t) :: centroids(this % n)
!
! For leaf return list of its points, for internal node return list of its
! centroids.
!
      integer :: i

      if (this % isleaf) then
        centroids = this % points(1:this % n)
      else
        do i=1, this % n
          centroids(i) = this % children(i) % p % centroid
        enddo
      endif 
    end function



!TODO make point/centroid the same type (point has radius component)??
    pure function ssnode_getradii(this) result(radii)
      class(ssnode_t), intent(in) :: this
      real(WP) :: radii(this % n)
!
! For internal node return list of its radii.
!
      integer :: i

      if (this % isleaf) then
        error stop 'getradii - node is a leaf'
      else
        do i=1, this % n
          radii(i) = this % children(i) % p % radius
        enddo
      endif 
    end function



    pure function mean(points, k)
      real(WP) :: mean
      integer, intent(in) :: k
      type(point_t), intent(in) :: points(:)

      if (size(points) < 1) &
          error stop 'mean - empty list of points'

      mean = sum(points % x(k)) / size(points)
    end function mean



    pure function variance(points, k)
      real(WP) :: variance
      integer, intent(in) :: k
      type(point_t), intent(in) :: points(:)

      real(WP) :: m

      m = mean(points, k)
      variance = sum((m - points % x(k))**2) / size(points)
    end function variance



    recursive subroutine validate(node, root, isvalid, level, points, nch)
      type(ssnode_t), pointer, intent(in) :: node
      type(ssnode_t), pointer, intent(in) :: root
      logical, intent(out) :: isvalid
      integer, intent(out) :: level, nch
      type(point_t), intent(out), allocatable :: points(:)

      logical :: isvalid0
      integer :: i, level0, nch0
      type(point_t), allocatable :: points0(:)

      if (.not. associated(node)) &
        error stop 'validate - unassocated node'

      ! no. of children must always be between M0 and M1
      isvalid = .false.
      if (node % n > M1) then
        print *, 'validation fails, too many children'
      elseif (node % n < M0 .and. .not. associated(node, root)) then
        print *, 'validation fails, too few children'
      elseif (node % n < 1) then
        print *, 'validation fails, too few children (even for a root)'
      else
        isvalid = .true.
      endif
      if (.not. isvalid) return

      ! leaf is valid if its points lie within its bounding envelope
      if (node % isleaf) then
        nch = 1
        level = 0
        points = node % getCentroids()
        do i = 1, size(points)
          if (points(i) % distance(node % centroid) <= node % radius) cycle 
          isvalid = .false.
          print *, 'validation fails, point is outside bounding envelope'
          exit
        enddo
        return
      endif

      ! internal node is valid if all its children are valid
      ! and on the same level...
      level = -1
      nch = 0
      allocate(points(0))
      do i = 1, node % n
        call validate(node%children(i)%p, root, isvalid0, level0, points0, &
        &             nch0)
        if (level == -1) level = level0
        if (.not. isvalid0) then
          isvalid = .false.
          print *, 'validation fails, one of childrens is not valid'
          return
        elseif (level /= level0) then
          isvalid = .false.
          print *, 'validation fails, childrens not at the same level'
          return
        endif
        points = [points, points0]
        nch = nch + nch0
      enddo
      level = level + 1
      nch = nch + 1

      ! ...and also if all its points lie within its bounding envelope
      do i = 1, size(points)
        if (points(i) % distance(node % centroid) <= node % radius) cycle 
        isvalid = .false.
        print *, 'validation fails, point is outside bounding envelope'
        exit
      enddo
    end subroutine validate



    logical function sstree_isvalid(this) result(isvalid)
      class(sstree_t), intent(in) :: this

      integer :: level, nch
      type(point_t), allocatable :: points(:)

      if (associated(this % root)) then
        call validate(this % root, this % root, isvalid, level, points, nch)
        print '(a,l1,a,i0,a,i0,a,i0)', &
            '(isvalid) valid? ',isvalid,'  level=',level,'  n=',size(points),&
            ' children=',nch
      else
        isvalid = .true.
        print '(a)', '(isvalid) empty tree'
      endif
    end function sstree_isvalid



    recursive subroutine finalize(node) 
      type(ssnode_t), pointer :: node

      integer :: i
      if (.not. associated(node)) then
        print *, 'finalize empty node'
        return
      endif

      ! recursively finalize children of internal nodes first...
      if (.not. node % isleaf) then
        do i=1, node % n
          call finalize(node % children(i) % p)
        enddo
      endif

      ! ...and deallocate the node itself
      deallocate(node)
    end subroutine finalize



    subroutine sstree_finalize(this)
      type(sstree_t), intent(inout) :: this
  print *, 'sstree_finalize...'
      call finalize(this % root)
  print *, '...sstree_finalize'
    end subroutine sstree_finalize

  end module sstree
