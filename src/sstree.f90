! Similarity Search (SS) Tree
! Ref: La Rocca, M., Advanced Algorithms and Data Structures, Ch10
!
! Tree constants should be M0 = k     k >= 2
!                          M1 = 2k-1
  module sstree
    use iso_fortran_env, only : output_unit
    implicit none
    private
    public sstree_t, point_t

    ! Dimensionality of ss-tree and stored points
    ! TODO temporarily made as module global parameters (to be moved to tree)
    !integer, parameter :: DIM=3, M0=2, M1=3
    integer, parameter :: DIM=3, M0=6, M1=6*2-1
    integer, parameter :: WP = kind(1.0d0)
    integer, parameter :: I8 = selected_int_kind(18)

    ! Debugging logs on/off
    !integer, parameter :: MSG_MEMORY = output_unit
    integer, parameter :: MSG_MEMORY = 0
    !integer, parameter :: MSG_DELETE = output_unit
    integer, parameter :: MSG_DELETE = 0

    ! Point class (public?)
    type point_t
      real(WP) :: x(DIM)
    contains
      procedure :: distance => point_distance
      procedure :: isin => point_isin
      procedure, private :: point_equals
      generic :: operator(==) => point_equals
    end type

    ! Array of ss-nodes
    type ssnode_ptr
      type(ssnode_t), pointer :: p => null()
    end type

    ! ss-node class (either internal node or a leaf)
    type ssnode_t
      type(point_t)    :: centroid
      real(WP)         :: radius
      type(ssnode_ptr) :: children(M1+1)
      type(point_t)    :: points(M1+1)
      logical          :: isleaf ! ... or internal node with children
      integer          :: n = 0
    contains
      procedure :: intersectsPoint => ssnode_intersectsPoint
      procedure :: findClosestChild => ssnode_findClosestChild
      procedure :: updateBoundingEnvelope => ssnode_updateBoundingEnvelope
      procedure :: directionOfMaxVariance => ssnode_directionOfMaxVariance
      procedure :: findSplitIndex => ssnode_findSplitIndex
      procedure :: minVarianceSplit => ssnode_minVarianceSplit
      procedure :: addPoint => ssnode_addPoint
      procedure :: deletePoint => ssnode_deletePoint
      procedure :: addChild => ssnode_addChild
      procedure :: deleteChild => ssnode_deleteChild
      procedure :: getCentroids => ssnode_getCentroids
      procedure :: getRadii => ssnode_getRadii
      procedure :: split => ssnode_split
      procedure :: siblingsToBorrowFrom => ssnode_siblingsToBorrowFrom
      procedure :: borrowFromSibling => ssnode_borrowFromSibling
      procedure :: findSiblingToMergeTo => ssnode_findSiblingToMergeTo
      procedure :: mergeChildren => ssnode_mergeChildren
    end type

    interface ssnode_t
      module procedure new_node
      module procedure new_leaf
    end interface

    ! ss-tree class (public)
    type sstree_t
      type(ssnode_t), pointer :: root => null()
      integer :: m0, m1, dim ! TODO not used at the moment
    contains
      procedure :: insert => sstree_insert
      procedure :: delete => sstree_delete
      procedure :: isvalid => sstree_isvalid
      procedure :: nearestNeighbor => sstree_nearestNeighbor
      final :: sstree_finalize
    end type

    interface sstree_t
      module procedure new_tree
    end interface

  contains

!
! Ss-node and ss-tree class constructors
!
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



!
! Ss-tree class finalizers
! 
    subroutine sstree_finalize(this)
      type(sstree_t), intent(inout) :: this
      if (MSG_MEMORY /= 0) &
        write(MSG_MEMORY,'(a)', advance='no') 'sstree_finalize...'
      call finalize(this % root)
      if (MSG_MEMORY /= 0) &
        write(MSG_MEMORY,'(a)') '   done.'
    end subroutine sstree_finalize

    recursive subroutine finalize(node) 
      type(ssnode_t), pointer :: node

      integer :: i
      if (.not. associated(node)) then
        print *, '(finalize) called with empty node'
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



!
! Point class methods
!
    elemental function point_distance(p0, p1) result(dis)
      real(WP) :: dis
      class(point_t), intent(in) :: p0, p1

      ! Euclidean distance between points P0 and P1
      dis = sum((p0 % x - p1 % x)**2)
      dis = sqrt(dis)
    end function point_distance



    elemental logical function point_equals(p0, p1)
      class(point_t), intent(in) :: p0, p1
      point_equals = all(p0 % x == p1 % x)
    end function point_equals



    pure integer function point_isin(this, arr)
      class(point_t), intent(in) :: this, arr(:)

      integer :: i

      point_isin = -1
      do i = 1, size(arr)
        if (arr(i) == this) then
          point_isin = i
          exit
        endif
      enddo
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



!
! Ss-nodes | points geometricaly-based helper functions
!

! Listing 10.10.
    elemental subroutine ssnode_updateBoundingEnvelope(this)
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



! Listing 10.4
    elemental logical function ssnode_intersectsPoint(this, point) result(res)
      class(ssnode_t), intent(in) :: this
      class(point_t), intent(in) :: point

      res = this % centroid % distance(point) <= this % radius
    end function ssnode_intersectsPoint



!
! Operations with ss-node containers (children / points arrays)
!
    subroutine ssnode_addPoint(this, point)
      class(ssnode_t), intent(inout) :: this
      class(point_t), intent(in) :: point

      if (.not. this % isleaf) &
          error stop 'addPoint - is not leaf'
      if (this % n > M1) &
          error stop 'addPoint - size oveflows'
      if (point % isin(this % points(1:this % n))>0) &
          error stop 'addPoint - duplicate point'

      this % n = this % n + 1
      this % points(this % n) = point
    end subroutine



    subroutine ssnode_deletePoint(this, point, ind)
      class(ssnode_t), intent(inout) :: this
      class(point_t), intent(in) :: point
      integer, intent(in), optional :: ind

      integer :: ind0, i

      if (.not. this % isleaf) &
          error stop 'deletePoint - is not leaf'
      if (this % n < 1) &
          error stop 'deletePoint - is empty'
      if (present(ind)) then
        ind0 = ind
      else
        ind0 = point % isin(this % points(1:this % n))
      endif
      if (ind0 < 1 .or. ind0 > this % n) &
          error stop 'deletePoint - point not present or wrong ind given'
      if (.not. this % points(ind0) == point) &
          error stop 'deletePoint - point and index are inconsistent'

      this % n = this % n - 1
      do i = ind0, this % n
        this % points(i) = this % points(i+1)
      enddo
    end subroutine



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



    subroutine ssnode_deleteChild(this, child, entry_only)
      class(ssnode_t), intent(inout) :: this
      class(ssnode_t), pointer, intent(inout) :: child
      logical, optional :: entry_only 
!
! entry_only == .false. (default): remove entry and deallocate node
! entry_only == .true.           : remove entry, leave node allocated
!
      integer :: i, i0
      logical :: entry_only0

      entry_only0 = .false. ! remove entry and deallocate on default
      if (present(entry_only)) entry_only0 = entry_only

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

      if (MSG_MEMORY /= 0 .and. .not. entry_only0) then
        write(MSG_MEMORY, '(a,i0,a,l1)') &
        '(deletechild) nold=', child%n, ' isleaf=', child%isleaf
      endif
      if (.not. entry_only0) deallocate(child)
    end subroutine



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
    !TODO make point/centroid the same type (point has radius component)??



    pure function sortIndices(keys) result(ind)
      real(WP), intent(in) :: keys(:)
      integer :: ind(size(keys))
      integer(I8) :: ibig
!
! Indirectly sort "key" array. Using insertion sort at the moment.
!
      integer :: i, j, n, i0
      real(WP) :: key0

      n = size(keys)
      ind = [(i, i=1,n)]

      do i = 2, n
        i0 = ind(i)
        key0 = keys(ind(i))
        j = i - 1
        do while ( j >= 1)
          if (keys(ind(j)) <= key0) exit
          ind(j+1) = ind(j)
          j = j-1
        enddo
        ind(j+1) = i0
      enddo

      !TODO temporarily verify 
      do i = 2, n
        if (keys(ind(i-1)) > keys(ind(i))) &
          error stop 'sortIndices - algorithm error'
      enddo 
      ibig = (n+n*n)/2 ! Gaussian: sum = (F+L)*n/2
      if (sum(int(ind,kind=I8)) /= ibig) error stop 'sum not correct' 

    end function sortIndices



!
! Validation method for ss-tree class
!
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
        error stop 'validate - unassociated node'
!print *, '(validate) isleaf=',node%isleaf,node%n

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
        if (.not. all(node % intersectsPoint(points))) then
          isvalid = .false.
          print *, 'validation fails, point is outside bounding envelope'
        endif

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
      if (.not. all(node % intersectsPoint(points))) then
        isvalid = .false.
        print *, 'validation fails, point is outside bounding envelope'
      endif
    end subroutine validate



!
! Insert method (and its helper functions)
!

!Listing 10.7
    subroutine sstree_insert(this, point)
      class(sstree_t) :: this
      type(point_t), intent(in) :: point
!
! Insert point to SS tree. Insert first point straight away, use recursive
! "insert" function otherwise.
! 
      type(ssnode_ptr) :: new(2)

      if (.not. associated(this % root)) then
        this % root => ssnode_t([point])
      else
        new = insert(this % root, point)
        if (associated(new(1) % p)) then
          if (MSG_MEMORY /= 0) &
          &   write(MSG_MEMORY, '(a,i0,a,l1)') &
          &   '(insert dealloc root) nold=', this%root%n, &
          &   ' isleaf=', this%root%isleaf
          deallocate(this % root)
          this % root => ssnode_t([new(1), new(2)])
        endif
      endif
    end subroutine sstree_insert



!Listing 10.6
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
        if (point % isin(node % points(1:node % n))>0) then
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
    subroutine ssnode_findSplitIndex(this, ind)
      class(ssnode_t), intent(inout) :: this
      integer, intent(out) :: ind
!
! Prepare for splitting the node. Reorder points/children and return
! index where to split.
!
      integer :: dir, i, j
      real(WP) :: key
      type(point_t) :: tmp_point
      type(ssnode_t), pointer :: tmp_node
      !integer, allocatable :: sind(:)
      !real(WP), allocatable :: keys(:)

      dir = this % directionOfMaxVariance()

      ! insertion sort for points or node pointers
      if (this % isleaf) then
        !sind = sortIndices(this % points(1:this % n) % x(dir))
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
        !allocate(keys(this % n))
        !do i = 1, this % n
        !  keys(i) = this % children(i) % p % centroid % x(dir)
        !enddo
        !sind = sortIndices(keys)
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

      ind = this % minVarianceSplit(dir)
    end subroutine ssnode_findSplitIndex



!Listing 10.9
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



!Listing 10.14
    pure function ssnode_minVarianceSplit(this, dir) result(split_index)
      integer :: split_index
      class(ssnode_t), intent(in) :: this
      integer, intent(in) :: dir
! 
! Get spliting index (index where the second half begins).
! Points/children arrays must be already sorted along dir's dimension
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




!Listing 10.5
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



 !Listing 10.8
    function ssnode_findClosestChild(this, target) result(res)
      class(ssnode_t), pointer :: res
      class(ssnode_t), intent(in) :: this
      type(point_t), intent(in) :: target
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
        dis = this % children(i) % p % centroid % distance(target)
        if (dis >= mindis) cycle
        mindis = dis
        res => this % children(i) % p
      enddo

      if (.not. associated(res)) & ! TODO can be removed after debug
          error stop 'findClosestChild found nothing'
    end function ssnode_findClosestChild



!
! Delete method (and its helper functions)
!
    subroutine sstree_delete(this, target)
      class(sstree_t), intent(inout) :: this
      type(point_t), intent(in) :: target

      logical :: info(2)
      integer, parameter :: DELETED=1, NEEDS_FIXING=2
      class(ssnode_t), pointer :: tmp_node

      if (.not. associated(this % root)) &
          error stop 'delete - empty tree'

      ! Recursive function to delete point and fix tree
      info = delete(this % root, target)

      ! Unsucessfull search
      if (.not. info(DELETED)) then
        print *, 'delete error - target not found'
        stop ! TODO deal with unsuccessfull search exception
      endif

      ! Does tree root pointer need updating?
      if (.not. info(NEEDS_FIXING)) return

      ! For level=0 graph, any number of points at root node is ok.
      ! Just free the node and nullify root after last point deleted.
      if (this % root % isleaf) then
        if (this % root % n ==0) then
          deallocate(this % root)
          if (MSG_DELETE /= 0) write(MSG_DELETE, '(a)') &
          & '(delete root) tree is now empty'
        endif

      else 
        ! If non-leaf root node has just one child, it can be freed and
        ! this child set as a new root (level of tree decreases)
        if (this % root % n == 1) then
          if (MSG_DELETE /= 0) write(MSG_DELETE, '(a)') &
          & '(delete root) tree level decreasing'
          tmp_node => this % root
          this % root => this % root % children(1) % p
          deallocate(tmp_node)
        elseif (.not.(this%root%n > 1 .and. this%root%n < M0)) then
          error stop 'delete tree - this branch should not be called'
        endif
      endif
    end subroutine sstree_delete



!Listing 10.15
    recursive function delete(node, target) result(info) 
      type(ssnode_t), pointer, intent(in) :: node
      type(point_t) :: target
      logical :: info(2) ! [was deleted?, needs fixing upstream?]

      integer, parameter :: DELETED=1, NEEDS_FIXING=2
      integer :: i
      class(ssnode_t), pointer :: nodeToFix
      type(ssnode_ptr), allocatable :: siblings(:)

      if (.not. associated(node)) &
        error stop 'delete - unassociated node'

      if (node % isleaf) then
        ! If leaf contains target, delete it. 
        ! Otherwise keep searching in the upstream instance.
        i = target % isin(node % points(1:node % n))  
        if (i > 0) then 
          call node % deletePoint(target, i)
          info(DELETED) = .true.        
          info(NEEDS_FIXING) = node % n < M0 
 !if (node % n > 0) call node % updateBoundingEnvelope()
 !TODO should we update envelope of leaf nodes?
        else 
          if (MSG_DELETE /= 0) write(MSG_DELETE, '(a)') &
          & '(delete) target not in leaf, continue searching' 
          info(DELETED) = .false.
          info(NEEDS_FIXING) = .false. 
        endif
        return
      endif

      if (MSG_DELETE /= 0) write(MSG_DELETE, '(a,i0,a)') &
        '(delete) non-leaf node, n=', node % n,'...'

      ! Node is a non-leaf node. Search recursively among its children until
      ! target is found or no more nodes remains to be traversed.
      ! Skip nodes whose boundary envelope does not contain target.
      nodeToFix => null()
      info(DELETED) = .false.
      do i = 1, node % n
        if (.not. node % children(i) % p % intersectsPoint(target)) cycle
        info = delete(node % children(i) % p, target)
        if (info(NEEDS_FIXING)) nodeToFix => node % children(i) % p
        if (info(DELETED)) exit
      enddo

      ! If no violation of the minimum number of entries in one of childrens
      ! has been made, we can return. If target was deleted, update bounding
      ! envelope.
      if (.not. associated(nodeToFix)) then
        if (info(DELETED)) call node % updateBoundingEnvelope()
        info(NEEDS_FIXING) = .false.  
        if (MSG_DELETE /= 0) write(MSG_DELETE, '(a,2l,a,i0)') &
        & '...(delete) nothing needs fixing (d/f) ',info,' n=',node % n
        return
      endif

      ! One of current children (nodeToFix) needs fixing
      siblings = node % siblingsToBorrowFrom(nodeToFix)
      if (size(siblings) > 0) then
        ! borrow item from one of the siblings
        if (MSG_DELETE /= 0) write(MSG_DELETE, '(a,i0)') &
        & '...(delete) borrow item, siblings available =', size(siblings)
        call nodeToFix % borrowFromSibling(siblings)
      elseif (node % n == 1) then
        !if (MSG_DELETE /= 0) write(MSG_DELETE, '(a)') &
        !& '...(delete) single child remains, this must be a root node' 
        !call node % mergeChildren(nodeToFix, null())
        error stop 'delete - node with one children not removed by root'
      elseif (node % n > 1) then 
        ! merge children
        if (MSG_DELETE /= 0) write(MSG_DELETE, '(a,i0)') &
        & '...(delete) merge children, siblings available =', node%n-1 
        call node % mergeChildren(nodeToFix, &
        &                         node % findSiblingToMergeTo(nodeToFix))
      else
        error stop 'delete - impossible branch called'
      endif
      call node % updateBoundingEnvelope()
      info(DELETED) = .true.
      info(NEEDS_FIXING) = node % n < M0
    end function delete



    function ssnode_siblingsToBorrowFrom(this, fixed) result(siblings)
      type(ssnode_ptr), allocatable :: siblings(:)
      class(ssnode_t), intent(in) :: this
      class(ssnode_t), pointer, intent(in) :: fixed
!
! On of this node's children (fixed) is one item short. Make a list of all
! fixed's sibling with more than "M0" items that can borrow one of its items
! to fixed.
!
      integer :: i, nfound
      logical :: fixed_contained
      type(ssnode_ptr), allocatable :: tmplist(:)

!print *, '(sibling to borrow from) ...'
      if (this % isleaf) &
          error stop 'siblingsToBorrowFrom - node is leaf'

      if (fixed % n >= M0) &
          error stop 'siblingsToBorrowFrom - fixed node has enough items'

      allocate(tmplist(M1))
      nfound = 0
      fixed_contained = .false.
      do i = 1, this % n
        if (associated(this % children(i) % p, fixed)) then
          fixed_contained = .true.
          cycle
        endif
        if (this % children(i) % p % n <= M0) cycle
        nfound = nfound + 1
        tmplist(nfound) % p => this % children(i) % p
      enddo
      if (.not. fixed_contained) &
          error stop 'siblingsToBorrowFrom - fixed is not among childrens'
      allocate(siblings(nfound))
      siblings = tmplist(1:nfound)
!print *, '... found possible siblings =',nfound
    end function ssnode_siblingsToBorrowFrom

    

! Listing 10.17
    subroutine ssnode_borrowFromSibling(this, siblings)
      class(ssnode_t), intent(inout) :: this
      type(ssnode_ptr) :: siblings(:)

      integer :: closest_i, i
      !type(ssnode_t), pointer :: closest_s
      type(point_t) :: closest_point
      class(ssnode_t), pointer :: closest_child, closest_s

      if (size(siblings) < 1) &
          error stop 'borrowFromSiblings - empty list'
      do i = 1, size(siblings)
        if (.not. associated(siblings(i)%p)) &
          error stop 'borrowFromSiblings - unassociated node'
        if (siblings(i) % p % n <= M0) &
          error stop 'borrowFromSiblings - too few items for a sibling'
      enddo

      call findClosestEntryInNodesList(siblings, this, closest_i, closest_s)

      if (closest_s % isleaf) then
        closest_point = closest_s % points(closest_i)
        call this % addPoint(closest_point)
        call closest_s % deletePoint(closest_point, closest_i)
      else
        closest_child => closest_s % children(closest_i) % p
        call this % addChild(closest_child)
        call closest_s% deleteChild(closest_child, .true.)
      endif
      call this % updateBoundingEnvelope()
      call closest_s % updateBoundingEnvelope()
    end subroutine ssnode_borrowFromSibling



! Listing 10.16
    subroutine findClosestEntryInNodesList(siblings, targetNode, ind, node)
      type(ssnode_ptr), intent(in) :: siblings(:)
      class(ssnode_t),  intent(in) :: targetNode
      integer, intent(out) :: ind
      class(ssnode_t), pointer, intent(out) :: node

      integer :: i, j, minj
      real(WP) :: dist, mindist, mindist_node

      node => null()
      mindist_node = huge(mindist_node)
      do i = 1, size(siblings)

        mindist = huge(mindist)
        if (siblings(i) % p % isleaf) then
          do j = 1, siblings(i) % p % n
            dist = targetNode % centroid % distance( &
                   siblings(i) % p % points(j))
            if (dist >= mindist) cycle
            minj = j
            mindist = dist
          enddo
        else
          do j = 1, siblings(i) % p % n
            dist = targetNode % centroid % distance( &
                   siblings(i) % p % children(j) % p % centroid)
            if (dist >= mindist) cycle
            minj = j
            mindist = dist
          enddo
        endif

        if (mindist >= mindist_node) cycle
        node => siblings(i) % p
        mindist_node = mindist
        ind = minj
      enddo
    end subroutine



    function ssnode_findSiblingToMergeTo(this, fixed) result(node)
      class(ssnode_t), intent(in) :: this
      class(ssnode_t), pointer, intent(in) :: fixed
      class(ssnode_t), pointer :: node
!
! Delete heuristics: which node to merge to?
! - for now find the closest sibling
! - TODO other choice would be the lower overlap
!
      integer :: i, mini
      real(WP) :: dist, mindist
      logical :: fixed_contained

      if (this % isleaf) &
          error stop 'findSiblingToMergeTo - this must be nonleaf'

      if (this % n < 2) &
          error stop 'findSiblingToMergeTo - no siblings'

      node => null()
      mindist = huge(mindist)
      fixed_contained = .false.
      do i = 1, this % n
        if (associated(fixed, this % children(i) % p)) then
          fixed_contained = .true.
          cycle
        endif
        dist = fixed % centroid % distance( &
               this % children(i) % p % centroid)
         if (dist >= mindist) cycle
         mindist = dist
         node => this % children(i) % p
      enddo

      if (.not. fixed_contained) &
          error stop 'findSiblingsToMergeTo - fixed not one of childrens'

      if (.not. associated(node)) &
          error stop 'findSiblingsToMergeTo - null siblings returning' 
    end function



!Listing 10.18
    subroutine ssnode_mergeChildren(this, childa, childb)
      class(ssnode_t), intent(inout) :: this
      class(ssnode_t), pointer :: childa, childb

      type(ssnode_t), pointer :: new

      if (associated(childb)) then
        new => merge2(childa, childb)
        call this % deleteChild(childa)
        call this % deleteChild(childb)
        call this % addChild(new)
      else
print *, '(mergeChildren) - nothing to merge, must be at root'
      endif
    end subroutine ssnode_mergeChildren



    function merge2(a, b) result(new)
      class(ssnode_t), intent(in) :: a, b
      type(ssnode_t), pointer :: new

      type(point_t), allocatable :: points(:)
      type(ssnode_ptr), allocatable :: children(:)

      if (a % isleaf .and. b % isleaf) then
        points = [a % points(1:a % n), b % points(1:b % n)]
        new => ssnode_t(points)
      elseif (.not. a % isleaf .and. .not. b % isleaf) then
        children = [a % children(1:a % n), b % children(1:b % n)]
        new => ssnode_t(children)
      else
          error stop 'merge2 - operands are not of the same type'
      endif
      if (new % n /= a % n + b % n) &
          error stop 'merge2 - sum of items is not same'
    end function merge2



!
! Search methods
!

    function sstree_nearestNeighbor(this, target) result(nn)
      class(sstree_t), intent(in) :: this
      class(point_t), intent(in) :: target
      type(point_t) :: nn

      real(WP) :: nndist
      integer :: nvisit

      nndist = huge(nndist)
      if (associated(this % root)) &
        call nearestNeighbor(this % root, target, nndist, nn, nvisit)

      if (nndist == huge(nndist)) then
        print *, 'nn not found'
      endif
  print *, '(NN search) found ', nn%distance(target), nvisit
  print *, nn%x

    end function



!Listing 10.19
    recursive subroutine nearestNeighbor(node, target, nndist, nn, nvisit)
      class(ssnode_t), intent(in) :: node
      class(point_t), intent(in) :: target
      real(WP), intent(inout) :: nndist
      type(point_t), intent(inout) :: nn
      integer, intent(out) :: nvisit

      integer :: i, nvisit0
      real(WP) :: dist
      type(point_t), allocatable :: points(:)
      integer, allocatable :: ind(:)
      real(WP), allocatable :: distarr(:)

      if (node % isleaf) then
        do i = 1, node % n
          dist = node % points(i) % distance(target)
          if (dist < nndist) then
            nndist = dist
            nn = node % points(i)
 print *, '(NN search) NN updated', nn%distance(target), nn%x
          end if
        enddo
        nvisit = 1

      else
        nvisit = 0
        points = node % getCentroids()
        distarr = target % distance(points)
        ind = sortIndices(distarr)
        do i = 1, node % n
          associate(child => node % children(ind(i)) % p)
  !         if (child % centroid % distance(target) < nndist) then
            if (child % centroid % distance(target)-child%radius < nndist) then
              call nearestNeighbor(child, target, nndist, nn, nvisit0)
              nvisit = nvisit + nvisit0
            endif
          end associate
        enddo
      endif
    end subroutine nearestNeighbor



!Listing 10.3
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

  end module sstree
