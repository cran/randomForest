c       Copyright (C) 2001  Leo Breiman and Adele Cutler

c       This program is free software; you can redistribute it and/or
c       modify it under the terms of the GNU General Public License
c       as published by the Free Software Foundation; either version 2
c       of the License, or (at your option) any later version.

c       This program is distributed in the hope that it will be useful,
c       but WITHOUT ANY WARRANTY; without even the implied warranty of
c       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c       GNU General Public License for more details.
c
c  Modified by Andy Liaw and Matt Wiener:
c  The main program is re-written as a C function to be called from R.
c  All calls to the uniform RNG is replaced with R's RNG.  Some subroutines
c  not called are excluded.  Variables and arrays declared as double as 
c  needed.  Unused variables are deleted.
c
c     SUBROUTINE MAKEA

      subroutine makea(x,mdim,nsample,cat,isort,v,a,b,mred)
      implicit double precision(a-h,o-z)
      double precision x(mdim,nsample),v(nsample)
      integer cat(mdim),isort(nsample),a(mdim,nsample),
     1     b(mdim,nsample)
      
c submakea constructs the mdim x nsample integer array a.  if there are
c less than 32,000 cases, this can be declared integer*2, otherwise
c integer*4. For each numerical variable with values x(m,n),
c n=1,...,nsample, the x-values are sorted from lowest to highest.
c Denote these by xs(m,n).  Then a(m,n) is the case number in which
c xs(m,n) occurs. The b matrix is also contructed here.
c if the mth variable is categorical, then a(m,n) is the category of the
c nth case number.

      do 10 mvar=1,mred

         if (cat(mvar).eq.1) then
            do 20 n=1, nsample
               v(n)=x(mvar,n)
               isort(n)=n
 20         continue
            call quicksort(v,isort,1,nsample,nsample)
            
            
            
c     this sorts the v(n) in ascending order. isort(n) is the case number 
c     of that v(n) nth from the lowest (assume the original case numbers 
c     are 1,2,...).  
            
            
            do 35 n=1,nsample-1
               n1=isort(n)
               n2=isort(n+1)
               a(mvar,n)=n1
               if(n.eq.1) b(mvar,n1)=1
               if (v(n).lt.v(n+1)) then
                  b(mvar,n2)=b(mvar,n1)+1
               else
                  b(mvar,n2)=b(mvar,n1)
               endif
 35         continue
            a(mvar,nsample)=isort(nsample)

         else

            do 40 ncat=1,nsample
               a(mvar,ncat)=nint(x(mvar,ncat))
 40         continue

         endif
 10   continue

      do 50 n=1,nsample
         isort(n)=n
 50   continue
      end
      
      
c     SUBROUTINE PREP
      
      subroutine prep(cl,nsample,nclass,ipi,pi,pid,nc,wtt)
      implicit double precision(a-h,o-z)      
      double precision pi(nclass),pid(nclass),wtt(nsample)
      integer nc(nclass),cl(nsample)

      call zerv(nc,nclass)
      do n=1,nsample
         nc(cl(n))=nc(cl(n))+1
      end do

      if (ipi.eq.0) then
         do 20 j=1,nclass
            pi(j)=dble(nc(j))/nsample
 20      continue
      endif
      
      sp=0
      do j=1,nclass
         sp=sp+pi(j)
      end do
      do j=1,nclass
         pi(j)=pi(j)/sp
      end do
      
      do 30 j=1,nclass
         if(nc(j).ge.1) then
            pid(j)=pi(j)*nsample/nc(j)
         else
            pid(j)=0
         end if
         do n=1,nsample
            wtt(n)=pid(cl(n))
         end do
 30   continue

      end
      
c     SUBROUTINE MODA
      
      subroutine moda(a,nuse,nsample,mdim,cat,maxcat,
     1     ncase,jin,ta)
      implicit double precision(a-h,o-z)
      integer a(mdim,nsample),cat(mdim),jin(nsample),
     1     ncase(nsample),ta(nsample)
      
      nuse=0
      do n=1,nsample
         if(jin(n).eq.1) nuse=nuse+1
      end do
      
      do m=1,mdim
         k=1
         nt=1
         if(cat(m).eq.1) then
            do n=1,nsample
               if(jin(a(m,k)).eq.1) then
                  a(m,nt)=a(m,k)
                  k=k+1
               else
                  do j=1,nsample-k
                     if(jin(a(m,k+j)).eq.1) then
                        a(m,nt)=a(m,k+j)
                        k=k+j+1
                        goto 28
                     end if
                  end do
               end if
 28            continue
               nt=nt+1
               if(nt.gt.nuse) goto 31
            end do
 31         continue
         end if
      end do
      
      if(maxcat.gt.1) then
         k=1
         nt=1
         do n=1,nsample
            if(jin(k).eq.1) then
               ncase(nt)=k
               k=k+1
            else
               do j=1,nsample-k
                  if(jin(k+j).eq.1) then
                     ncase(nt)=k+j
                     k=k+j+1
                     goto 58
                  end if
               end do
            end if
 58         continue
            nt=nt+1
            if(nt.gt.nuse) goto 61
         end do
 61      continue
      end if
      end
      
      

      
c     SUBROUTINE BUILDTREE
      
      subroutine buildtree(a,b,cl,cat,mdim,nsample,nclass,treemap,
     1     bestvar,bestsplit,bestsplitnext,tgini, nodestatus,nodepop,
     1     nodestart,classpop,tclasspop,tclasscat,ta,nrnodes,
     1     idmove,ndsize,ncase,parent,jin,mtry,iv,nodeclass,ndbigtree,
     1     win,wr,wc,wl,mred,nuse,mind)
      

c Buildtree consists of repeated calls to two subroutines, Findbestsplit
c and Movedata.  Findbestsplit does just that--it finds the best split of
c the current node.  Movedata moves the data in the split node right and
c left so that the data corresponding to each child node is contiguous.
c The buildtree bookkeeping is different from that in Friedman's original
c CART program.  ncur is the total number of nodes to date.
c nodestatus(k)=1 if the kth node has been split.  nodestatus(k)=2 if the
c node exists but has not yet been split, and =-1 of the node is terminal.
c A node is terminal if its size is below a threshold value, or if it is
c all one class, or if all the x-values are equal.  If the current node k
c is split, then its children are numbered ncur+1 (left), and
c ncur+2(right), ncur increases to ncur+2 and the next node to be split is
c numbered k+1.  When no more nodes can be split, buildtree returns to the
c main program.
      
      implicit double precision(a-h,o-z)
      integer a(mdim,nsample),cl(nsample),cat(mdim),
     1     treemap(2,nrnodes),bestvar(nrnodes),
     1     bestsplit(nrnodes), nodestatus(nrnodes),ta(nsample),
     1     nodepop(nrnodes),nodestart(nrnodes),
     1     bestsplitnext(nrnodes),idmove(nsample),
     1     ncase(nsample),parent(nrnodes),b(mdim,nsample),
     1     jin(nsample),iv(mred),nodeclass(nrnodes),mind(mred)
      
      
      double precision tclasspop(nclass),classpop(nclass,nrnodes),
     1     tclasscat(nclass,32),win(nsample),wr(nclass),wc(nclass),
     1     wl(nclass),tgini(mdim)
      integer msplit
      
      msplit = 0
      call zerv(nodestatus,nrnodes)
      call zerv(nodestart,nrnodes)
      call zerv(nodepop,nrnodes)
      call zermr(classpop,nclass,nrnodes)
      
      
      
      do 20 j=1,nclass
         classpop(j,1)=tclasspop(j)
 20   continue

      ncur=1
      nodestart(1)=1
      nodepop(1)=nuse
      nodestatus(1)=2
      
c     start main loop

      do 30 kbuild=1,nrnodes

         if (kbuild.gt.ncur) goto 50
         if (nodestatus(kbuild).ne.2) goto 30
         
c     initialize for next call to findbestsplit

         ndstart=nodestart(kbuild)
         ndend=ndstart+nodepop(kbuild)-1
         do 40 j=1,nclass
            tclasspop(j)=classpop(j,kbuild)
 40      continue
         jstat=0

         call findbestsplit(a,b,cl,mdim,nsample,nclass,cat,ndstart,
     1        ndend,tclasspop,tclasscat,msplit,decsplit,nbest,ncase,
     1        jstat,jin,mtry,iv,win,wr,wc,wl,mred,kbuild,mind)
         
         
         if(jstat.eq.1) then
            nodestatus(kbuild)=-1
            goto 30
         else
            bestvar(kbuild)=msplit
            tgini(msplit)=decsplit+tgini(msplit)
            if (cat(msplit).eq.1) then
               bestsplit(kbuild)=a(msplit,nbest)
               bestsplitnext(kbuild)=a(msplit,nbest+1)
            else
               bestsplit(kbuild)=nbest
               bestsplitnext(kbuild)=0
            endif
         endif
                  
         call movedata(a,ta,mdim,nsample,ndstart,ndend,idmove,ncase,
     1        msplit,cat,nbest,ndendl)
                  
c     leftnode no.= ncur+1, rightnode no. = ncur+2.

         nodepop(ncur+1)=ndendl-ndstart+1
         nodepop(ncur+2)=ndend-ndendl
         nodestart(ncur+1)=ndstart
         nodestart(ncur+2)=ndendl+1

c     find class populations in both nodes
         
         do 60 n=ndstart,ndendl
            if(cat(msplit).gt.1)then
               nc=ncase(n)
            else
               nc=ncase(n)
            end if
            j=cl(nc)
            classpop(j,ncur+1)=classpop(j,ncur+1)+win(nc)
 60      continue
         do 70 n=ndendl+1,ndend
            if(cat(msplit).gt.1) then
               nc=ncase(n)
            else
               nc=ncase(n)
            end if
            j=cl(nc)
            classpop(j,ncur+2)=classpop(j,ncur+2)+win(nc)
 70      continue


c     check on nodestatus

         nodestatus(ncur+1)=2
         nodestatus(ncur+2)=2
         if (nodepop(ncur+1).le.ndsize) nodestatus(ncur+1)=-1
         if (nodepop(ncur+2).le.ndsize) nodestatus(ncur+2)=-1
         popt1=0
         popt2=0
         do j=1,nclass
            popt1=popt1+classpop(j,ncur+1)
            popt2=popt2+classpop(j,ncur+2)
         end do
         
         do j=1,nclass
            if (classpop(j,ncur+1).eq.popt1) nodestatus(ncur+1)=-1
            if (classpop(j,ncur+2).eq.popt2) nodestatus(ncur+2)=-1
         end do

         treemap(1,kbuild)=ncur+1
         treemap(2,kbuild)=ncur+2
         parent(ncur+1)=kbuild
         parent(ncur+2)=kbuild
         nodestatus(kbuild)=1
         ncur=ncur+2
         if (ncur.ge.nrnodes) goto 50
         
 30   continue
 50   continue

      ndbigtree=nrnodes
      do k=nrnodes,1,-1
         if (nodestatus(k).eq.0) ndbigtree=ndbigtree-1
         if (nodestatus(k).eq.2) nodestatus(k)=-1
      end do

      
      do kn=1,ndbigtree
         if(nodestatus(kn).eq.-1) then
            pp=0
            do j=1,nclass
               if(classpop(j,kn).gt.pp) then
                  nodeclass(kn)=j
                  pp=classpop(j,kn)
               end if
            end do
         end if
      end do
            
      end


c     SUBROUTINE FINDBESTSPLIT

c For the best split, msplit is the variable split on. decsplit is the
c dec. in impurity.  If msplit is numerical, nsplit is the case number
c of value of msplit split on, and nsplitnext is the case number of the
c next larger value of msplit.  If msplit is categorical, then nsplit is
c the coding into an integer of the categories going left.

      subroutine findbestsplit(a,b,cl,mdim,nsample,nclass,cat,
     1     ndstart,ndend,tclasspop,tclasscat,msplit,decsplit,nbest,
     1     ncase,jstat,jin,mtry,iv,win,wr,wc,wl,mred,kbuild,mind)
      implicit double precision(a-h,o-z)      
      integer a(mdim,nsample),cl(nsample),cat(mdim),iv(mred),
     1     ncase(nsample),b(mdim,nsample),jin(nsample), nn, j          
      double precision tclasspop(nclass),tclasscat(nclass,32),
     1     win(nsample),
     1     wr(nclass),wc(nclass),wl(nclass), xrand
      integer mind(mred)
c      real pno, pdo, rrd, rld
      
c     compute initial values of numerator and denominator of Gini
      
      pno=0.0
      pdo=0.0
      do 10 j=1,nclass
         pno=pno+tclasspop(j)*tclasspop(j)
         pdo=pdo+tclasspop(j)
 10   continue
      crit0=pno/pdo
      jstat=0
c      zz=rrand()
c      call rrand(zz)
            
c     start main loop through variables to find best split
      
      critmax=-1.0e20
      
      do k = 1, mred
         mind(k) = k
      end do

      nn = mred
      do 20 mt=1,mtry
 200     continue
c
c  sampling mtry variables w/o replacement.
c
         call rrand(xrand)
         j = int(nn * xrand) + 1
         mvar = mind(j)
         mind(j) = mind(nn)
         nn = nn - 1
c
         if(cat(mvar).eq.1) then
            rrn=pno
            rrd=pdo
            rln=0
            rld=0
            call zervr(wl,nclass)
            do 50 j=1,nclass
               wr(j)=tclasspop(j)
 50         continue
            critvar=-1e20
            
            do 60 nsp=ndstart,ndend-1
               nc=a(mvar,nsp)
               u=win(nc)
               k=cl(nc)
               rln=rln+u*(2*wl(k)+u)
               rrn=rrn+u*(-2*wr(k)+u)
               rld=rld+u
               rrd=rrd-u
               wl(k)=wl(k)+u
               wr(k)=wr(k)-u
                              
               if (b(mvar,nc).lt.b(mvar,a(mvar,nsp+1))) then
                  if(dmin1(rrd,rld).gt.1.0e-5) then
                     crit=(rln/rld)+(rrn/rrd)
                     if (crit.gt.critvar) then
                        nbestvar=nsp
                        critvar=crit
                     endif
                  end if
               end if
 60         continue
 65         continue

            if (critvar.gt.critmax) then
               msplit=mvar
               nbest=nbestvar
               critmax=critvar
            endif
            
         else

c     compute the decrease in impurity given by categorical splits

            lcat=cat(mvar)
            call zermr(tclasscat,nclass,32)
            do 70 nsp=ndstart,ndend
               nc=ncase(nsp)
               l=a(mvar,ncase(nsp))
               tclasscat(cl(nc),l)=tclasscat(cl(nc),l)+win(nc)
 70         continue
            nnz=0
            do i=1,lcat
               su=0
               do j=1,nclass
                  su=su+tclasscat(j,i)
               end do
               if(su.gt.0) nnz=nnz+1
            end do
            if (nnz.eq.1) then
               critvar=-1.0e25
            else
               call catmax(pno,pdo,tclasscat,tclasspop,nclass,lcat,
     1              nbestvar,critvar)
            end if
                        
c this last subroutine returns those categories going left in the best split. 
c This is coded into a long integer (see under subroutine catmax below for 
c details). 
          
            if (critvar.gt.critmax) then
               msplit=mvar
               nbest=nbestvar
               critmax=critvar
            endif
         endif
 20   continue
 25   continue                
      decsplit=critmax-crit0
      if (critmax.lt.-1.0e10) jstat=1
      
      end
      
C     SUBROUTINE CATMAX

      subroutine catmax(pno,pdo,tclasscat,tclasspop,nclass,lcat,
     1     ncatsplit,rmaxdec)
            
c this subroutine finds the best categorical split of a categorical variable
c with lcat categories, nclass classes and tclasscat(j,l) is the number of 
c cases in class j with category value l. The method used is an exhaustive 
c search over all partitions of the category values.  For the two class 
c problem, there is a faster exact algorithm we will add later.  
c If lcat.ge.10, the exhaustive search gets slow and there is a faster 
c iterative algorithm we can add later.
      
      parameter(jmax=100)
      implicit double precision(a-h,o-z)
      double precision tclasscat(nclass,32),tclasspop(nclass),
     1     tmpclass(jmax)
      integer icat(32) 
      
      rmaxdec=-1e20
      do 10 n=1,(2**(lcat-1))-1
         call myunpack(lcat,n,icat)
         call zervr(tmpclass,jmax)
         do 20 l=1,lcat
            if(icat(l).eq.1) then
               do 30 j=1, nclass
                  tmpclass(j)=tmpclass(j)+tclasscat(j,l)
 30            continue
            endif
 20      continue
         pln=0
         pld=0
         do 40 j=1,nclass
            pln=pln+tmpclass(j)*tmpclass(j)
            pld=pld+tmpclass(j)
 40      continue
         prn=0
         do 50 j=1,nclass
            tmpclass(j)=tclasspop(j)-tmpclass(j)
            prn=prn+tmpclass(j)*tmpclass(j)
 50      continue
         tdec=(pln/pld)+((prn)/(pdo-pld))
         if (tdec.gt.rmaxdec) then
            rmaxdec=tdec
            ncatsplit=n
         endif
 10   continue
      
      end

      
c     SUBROUTINE MOVEDATA     

c This subroutine is the heart of the buildtree construction. Based on the 
c best split the data in the part of the a matrix corresponding to the 
c current node is moved to the left if it belongs to the left child and 
c right if it belongs to the right child.

      subroutine movedata(a,ta,mdim,nsample,ndstart,ndend,idmove,
     1     ncase,msplit,cat,nbest,ndendl)
      implicit double precision(a-h,o-z)
      integer a(mdim,nsample),ta(nsample),idmove(nsample),
     1     ncase(nsample),cat(mdim),icat(32)
      
c     compute idmove=indicator of case nos. going left
      
      if (cat(msplit).eq.1) then
         do 10 nsp=ndstart,nbest
            nc=a(msplit,nsp)
            idmove(nc)=1
 10      continue
         do 20 nsp=nbest+1, ndend
            nc=a(msplit,nsp)
            idmove(nc)=0
 20      continue
         ndendl=nbest
      else
         ndendl=ndstart-1
         l=cat(msplit)
         call myunpack(l,nbest,icat)
         do 30 nsp=ndstart,ndend
            nc=ncase(nsp)
            if (icat(a(msplit,nc)).eq.1) then
               idmove(nc)=1
               ndendl=ndendl+1
            else
               idmove(nc)=0
            endif
 30      continue
      endif
      
c     shift case. nos. right and left for numerical variables.        
      
      do 40 msh=1,mdim
         if (cat(msh).eq.1) then
            k=ndstart-1
            do 50 n=ndstart,ndend
               ih=a(msh,n)
               if (idmove(ih).eq.1) then
                  k=k+1
                  ta(k)=a(msh,n)
               endif
 50         continue
            do 60 n=ndstart,ndend
               ih=a(msh,n)
               if (idmove(ih).eq.0) then 
                  k=k+1
                  ta(k)=a(msh,n)
               endif
 60         continue

            do 70 k=ndstart,ndend
               a(msh,k)=ta(k)
 70         continue
         endif

 40   continue
      ndo=0
      if(ndo.eq.1) then
	 do 140 msh=1,mdim
            if (cat(msh).gt.1) then
               k=ndstart-1
               do 150 n=ndstart,ndend
                  ih=ncase(n)
                  if (idmove(ih).eq.1) then
                     k=k+1
                     ta(k)=a(msh,ih)
                  endif
 150           continue
               do 160 n=ndstart,ndend
                  ih=ncase(n)
                  if (idmove(ih).eq.0) then 
                     k=k+1
                     ta(k)=a(msh,ih)
                  endif
 160           continue

               do 170 k=ndstart,ndend
                  a(msh,k)=ta(k)
 170           continue
            endif

 140     continue
      end if

c     compute case nos. for right and left nodes.
      
      if (cat(msplit).eq.1) then
         do 80 n=ndstart,ndend
            ncase(n)=a(msplit,n)
 80      continue
      else
         k=ndstart-1
         do 90 n=ndstart, ndend
            if (idmove(ncase(n)).eq.1) then
               k=k+1
               ta(k)=ncase(n)
            endif
 90      continue
         do 100 n=ndstart,ndend
            if (idmove(ncase(n)).eq.0) then
               k=k+1
               ta(k)=ncase(n)
            endif
 100     continue
         do 110 k=ndstart,ndend
            ncase(k)=ta(k)
 110     continue
      endif
      
      end
      

c     SUBROUTINE XTRANSLATE

c this subroutine takes the splits on numerical variables and translates them
c back into x-values.  It also unpacks each categorical split into a 
c 32-dimensional vector with components of zero or one--a one indicates that 
c the corresponding category goes left in the split.

      subroutine xtranslate(x,mdim,nrnodes,nsample,bestvar,
     1     bestsplit,bestsplitnext,xbestsplit,nodestatus,cat,
     1     ndbigtree)
      implicit double precision(a-h,o)      
      integer cat(mdim),bestvar(nrnodes),bestsplitnext(nrnodes),
     1     nodestatus(nrnodes),bestsplit(nrnodes)
      double precision x(mdim,nsample),xbestsplit(nrnodes)

      do 10 k=1,ndbigtree
         if (nodestatus(k).eq.1) then
            m=bestvar(k)
            if (cat(m).eq.1) then
               xbestsplit(k)=(x(m,bestsplit(k))+
     1              x(m,bestsplitnext(k)))/2
            else
               xbestsplit(k)=dble(bestsplit(k))
            endif
         endif
 10   continue
      end

C     SUBROUTINE TESTREEBAG
      

      subroutine testreebag(xts,nts,mdim,treemap,nodestatus,
     1     xbestsplit,cbestsplit,bestvar,nodeclass,nrnodes,
     1     ndbigtree,cat,nclass,jts,nodex,maxcat)
      
      implicit double precision (a-h,o-z)
      double precision xts(mdim,nts),xbestsplit(nrnodes)
      
      integer treemap(2,nrnodes),bestvar(nrnodes),
     1     nodeclass(nrnodes),cat(mdim),nodestatus(nrnodes),jts(nts),
     1     nodex(nts),cbestsplit(maxcat,nrnodes),icat(32)
      integer k, l, ncat, j, n, kt, jcat
      
      call zerv(jts,nts)
      call zerv(nodex,nts)
      
      do k=1,ndbigtree
         if(bestvar(k).gt.0) l=cat(bestvar(k))
         if(l.gt.1) then
            ncat=nint(xbestsplit(k))
            call myunpack(l,ncat,icat)
            do j=1,l
               cbestsplit(j,k)=icat(j)
            end do
         end if
      end do
            
      do n=1,nts
         kt=1
         do k=1,ndbigtree
            if (nodestatus(kt).eq.-1) then
               jts(n)=nodeclass(kt)
               nodex(n)=kt
               goto 100
            end if
            m=bestvar(kt)
            if (cat(m).eq.1) then
               if (xts(m,n).le.xbestsplit(kt)) then 
                  kt=treemap(1,kt)
               else
                  kt=treemap(2,kt)
               endif
               
            else
               jcat=nint(xts(m,n))
               if (cbestsplit(jcat,kt).eq.1) then
                  kt=treemap(1,kt)
               else
                  kt=treemap(2,kt)
               endif
            endif
         end do
 100     continue
      end do
      
      end
      
C     SUBROUTINE COMPTSERR
      
      subroutine comptserr(countts,jts,clts,jet,ntest,nclass,errts,
     1     pid,labelts)
      implicit double precision (a-h,o-z)
      integer jts(ntest),clts(ntest),jet(ntest), cmax, ntest, nclass
      double precision countts(nclass,ntest), pid(nclass)
      
      rmissts=0.0
      do n=1,ntest
         countts(jts(n),n)=countts(jts(n),n)+1
      end do

      do n=1,ntest
         cmax=0
         do j=1,nclass
            if (countts(j,n).gt.cmax) then
               jet(n)=j
               cmax=countts(j,n)
            end if
         end do
      end do

      if(labelts.eq.1) then
         do n=1,ntest
            if (jet(n).ne.clts(n)) rmissts = rmissts + 1.0
         end do
         errts=dble(rmissts)/dble(ntest)
      end if
      end
      
C     SUBROUTINE UNIF & FUNCTION NSELECT
      
      subroutine createclass(x,cl,ns,nsample,mdim,tx,p,
     1     sm,ndble,iaddcl)
      implicit double precision (a-h,o-z)
      double precision x(mdim,nsample),tx(ns),sm(ns),p(ns), xrand
      integer ndble(ns),cl(nsample)
      
      do n=1,ns
         cl(n)=1
      end do
      do n=ns+1,nsample
         cl(n)=2
      end do
            
      if(iaddcl.eq.1) then
         do n=ns+1,nsample
            do m=1,mdim
               call rrand(xrand)
               k=int(xrand*ns)+1
               x(m,n)=x(m,k)
            end do
         end do
      end if
      
      if(iaddcl.eq.2) then
         do m=1,mdim
            do n=1,ns
               tx(n)=x(m,n)
            end do
            call quicksort(tx,ndble,1,ns,ns)
            do k=2,ns-1
               p(k)=.5*(tx(k+1)-tx(k-1))/(tx(ns)-tx(1))
            end do
            p(ns)=.5*(tx(2)-tx(1))/(tx(ns)-tx(1))
            
            sm(1)=p(1)
            do k=2,ns
               sm(k)=sm(k-1)+p(k)
            end do
            
            do n=ns+1,nsample
               k=nselect(ns,sm)
               x(m,n)=tx(k)
            end do
         end do 
      end if 
      end
      
      integer function nselect(ns,sm)
      double precision sm(ns), u
      nselect=0
      call rrand(u)
      kp=ns
      km=0
      do j=1,1000
         kt=nint(dble(kp+km)/2)
         if (u.lt.sm(kt)) then
            kp=kt
         else
            km=kt
         end if
         if (kp-km.eq.1) then
            nselect=kp
            goto 69
         end if
      end do
 69   return
      end
      
      
C     SUBROUTINE OOB
      
      subroutine oob(nsample,nclass,jin,cl,jtr,jerr,counttr,out,
     1     errtr,errc,rmargin,q,jest,wtt)
      implicit double precision (a-h,o-z)      
      integer jin(nsample),cl(nsample),jtr(nsample),out(nsample),
     1     jerr(nsample),jest(nsample),counttr(nclass,nsample)
      double precision rmargin(nsample),q(nclass,nsample),wtt(nsample),
     1     errtr, rmiss, rmissc
      
      inderr=1
      if(inderr.eq.1) then
         outc=0
         rmissc=0.0
         do n=1,nsample
            if(jin(n).eq.0) then
               outc=outc+1.
               if(jtr(n).ne.cl(n)) rmissc=rmissc+1.
            end if
         end do
         errc=100*rmissc/outc
      end if
      
      call zerv(jerr,nsample)
      
      rmiss=0.0
      do n=1,nsample
         if(out(n).gt.0) then
            smax=0.0
            smaxtr=0.0
            do j=1,nclass
               q(j,n)=dble(counttr(j,n))/out(n)
               if(j.ne.cl(n)) smax=dmax1(q(j,n),dble(smax))
               if (q(j,n).gt.smaxtr) then
                  smaxtr=q(j,n)
                  jest(n)=j
               end if
            end do
            if(jest(n).ne.cl(n)) then
               rmiss=rmiss+1.0
               jerr(n)=1
            end if
            pth=q(cl(n),n)
            rmargin(n)=pth-smax
         end if
      end do
      errtr=rmiss/nsample
      end
      
C     SUBROUTINE PERMOBAR
      
      subroutine permobmr(mr,x,tp,tx,jin,nsample,mdim)
      implicit double precision (a-h,o-z)
      dimension jin(nsample)
      double precision tp(nsample), x(mdim, nsample), tx(nsample)
      kout=0
      call zervr(tp,nsample)
      do n=1,nsample
         if(jin(n).eq.0) then
            kout=kout+1
            tp(kout)=x(mr,n)
         end if
      end do 
      call perm1(kout,nsample,tp)
      iout=0
      do n=1,nsample
         tx(n)=x(mr,n)
         if(jin(n).eq.0) then
            iout=iout+1
            x(mr,n)=tp(iout)
         end if
      end do
      end
      

      
C     SUBROUTINE FINISHIMP
C     Modified by A. Liaw 2/11/2002 (removed graph from argument)
      subroutine finishimp(rmissimp,countimp,out,cl,nclass,mdim,
     1     nsample, errimp,rimpmarg,diffmarg,cntmarg,
     1     rmargin,counttr,jest,errtr)
      implicit double precision (a-h,o-z)      
      integer cl(nsample),countimp(nclass,nsample,mdim),
     1     counttr(nclass,nsample),out(nsample),jest(nsample)
      
      double precision errimp(mdim),rimpmarg(mdim,nsample),
     1     rmissimp(mdim),errtr,
     1     diffmarg(mdim),cntmarg(mdim),rmargin(nsample)
c     1	graph(nclass,nsample,mdim)

      imax=0
      ks0=0
      call zervr(rmissimp,mdim)
      call zervr(errimp,mdim)
      
      do m1=1,mdim
         do n=1,nsample
            if(out(n).ge.1) then
               lmax=0
               kmax=0
               do j=1,nclass
                  ks=countimp(j,n,m1)+counttr(j,n)
                  if (ks.gt.kmax) then
                     kmax=ks
                     imax=j
                  end if
                  if(j.ne.cl(n)) lmax=max0(lmax,ks)
                  if(j.eq.cl(n)) ks0=ks
               end do
               if(imax.ne.cl(n)) rmissimp(m1)=rmissimp(m1)+1
               rimpmarg(m1,n)=dble(ks0-lmax)/out(n)
            end if
         end do 
      end do 
      do m1=1,mdim
         errimp(m1)=rmissimp(m1)/nsample
         errimp(m1)=100*(errimp(m1)-errtr)/errtr
         errimp(m1)=dmax1(0.0d0,errimp(m1))
      end do
      do m=1,mdim
         diffmarg(m)=0
         cntmarg(m)=0
         do n=1,nsample
            diffmarg(m)=diffmarg(m)+(rmargin(n)-rimpmarg(m,n))
            if(rimpmarg(m,n).lt.rmargin(n)) cntmarg(m)=cntmarg(m)+1 
            if(rimpmarg(m,n).gt.rmargin(n)) cntmarg(m)=cntmarg(m)-1
         end do
         diffmarg(m)=100*diffmarg(m)/nsample
         cntmarg(m)=cntmarg(m)/nsample
      end do
            
      end
      
C     SUBROUTINE LOCATEOUT
      
      subroutine locateout(prox,cl,near,nsample,nclass,ncp,
     1     iaddcl,outlier,tout,isort,clp)
      implicit double precision (a-h,o-z)      
      double precision prox(near,near)
      double precision outlier(near),tout(near)
      integer ncp(near),cl(nsample),isort(nsample),ntt(0:30),
     1     clp(near)
      
      call zervr(outlier,near)
      
      if (iaddcl.eq.1) then
         jpclass=1
      else
         jpclass=nclass
      end if
      ntt(0)=0
      nt=0
      do jp=1,jpclass
         
         do n=1,near
            if(cl(n).eq.jp) then
               nt=nt+1
               ncp(nt)=n
            end if
         end do
         ntt(jp)=nt
         n1=ntt(jp-1)+1
         n2=ntt(jp)
         
         do i=n1,n2
            outlier(i)=0
            do j=n1,n2
               if(j.ne.i) then
                  outlier(i)=outlier(i)+(dble(prox(ncp(i),ncp(j))))**2
               end if
            end do
         end do
         
         do i=n1,n2
            outlier(i)=1.0/outlier(i)
            tout(i)=outlier(i)
            clp(i)=jp
         end do
         
         call quicksort(tout,isort,n1,n2,nsample)
         rmed=tout((n1+n2)/2)
         dev=0
         do i=n1,n2
            dev=dev+abs(tout(i)-rmed)
         end do
         dev=dev/(n2-n1+1)
         do i=n1,n2
            outlier(i)=(outlier(i)-rmed)/dev
            outlier(i)=dmax1(outlier(i),0.0d0)
         end do
         
      end do
      
      end
      
C     SUBROUTINE RUNFOREST

c     reads a forest file and runs dnew data through it
c mdim is the number of predictors
c ntest is the number of data points to predict for
c clts is whether test set has labels
c xts is an array of doubles:  mdim x ntest
c countts is a blank array of doubles:  nclass x ntest
c cbestsplit  is a blank array of integers:  maxcat x nrnodes
c jts is an integer vector of length ntest
c jet is an integer vector of length ntest
c nodexts is an integer vector of length ntest
c prox is 0 if proximity is not wanted, 1 if it is
c proxmatrix is a matrix of dimension ntest x ntest
c   if prox is 1, and is unused (and a single double) otherwise
c
      subroutine runforest(mdim,ntest,nclass,maxcat,nrnodes,
     1     labelts,jbt,clts,xts,xbestsplit,pid,countts,treemap,
     1     nodestatus,cat,cbestsplit,nodeclass,jts,jet,bestvar,
     1     nodexts,ndbigtree, prox, proxmatrix)
      
      implicit double precision (a-h,o-z)
      double precision xts(mdim,ntest),xbestsplit(nrnodes,jbt),
     1     pid(nclass),
     1     countts(nclass,ntest),errts, proxmatrix(ntest, ntest)
      
      integer treemap(2,nrnodes,jbt),nodestatus(nrnodes,jbt),
     1     cat(mdim),cbestsplit(maxcat,nrnodes),
     1     nodeclass(nrnodes,jbt),
     1     bestvar(nrnodes,jbt),jts(ntest),clts(ntest),jet(ntest),
     1     nodexts(ntest),ndbigtree(jbt),prox,n1,n2

      call zermr(countts,nclass,ntest)
      do jb=1,jbt
         call testreebag(xts,ntest,mdim,treemap(1,1,jb),
     1        nodestatus(1,jb),
     1        xbestsplit(1,jb),cbestsplit,bestvar(1,jb),
     1        nodeclass(1,jb), 
     1        nrnodes, ndbigtree(jb),cat,nclass,jts,nodexts,maxcat)
c     
c if desired, do proximities for this round
c
         if (prox.eq.1) then
            do n1=1,ntest
               do n2=1,ntest
                  if(nodexts(n1).eq.nodexts(n2)) then
                     proxmatrix(n1,n2)=proxmatrix(n1,n2) + 1.0
                  end if
               end do
            end do
         end if
         call comptserr(countts,jts,clts,jet,ntest,nclass,
     1        errts,pid,labelts)
         
      end do
      
c     
c     if proximities requested, do the final adjustment
c     (division by 2 * number of trees)
c     
      
      if (prox.eq.1) then     
         do n1=1,ntest
            do n2=(n1+1),ntest
               proxmatrix(n1,n2)=proxmatrix(n1,n2)/jbt
               proxmatrix(n2,n1)=proxmatrix(n1,n2)
            end do
            proxmatrix(n1,n1) = 1.0
         end do
      end if
      return
      end




C     SUBROUTINE QUICKSORT
      
      subroutine quicksort (v,iperm,ii,jj,kk)
c     
c     puts into iperm the permutation vector which sorts v into
c     increasing order.  only elementest from ii to jj are considered.
c     array iu(k) and array il(k) permit sorting up to 2**(k+1)-1 elements
c     
c     this is a modification of acm algorithm #347 by r. c. singleton,
c     which is a modified hoare quicksort.
c     
      implicit double precision (a-h,o-z)
      dimension iperm(kk),iu(32),il(32)
      integer t,tt
      double precision v(kk)
      
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=iperm(ij)
      vt=v(ij)
      if (v(i).le.vt) go to 30
      iperm(ij)=iperm(i)
      iperm(i)=t
      t=iperm(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      iperm(ij)=iperm(j)
      iperm(j)=t
      t=iperm(ij)
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      iperm(ij)=iperm(i)
      iperm(i)=t
      t=iperm(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   iperm(l)=iperm(k)
      iperm(k)=tt
      v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      tt=iperm(l)
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=iperm(i+1)
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  iperm(k+1)=iperm(k)
      v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      iperm(k+1)=t
      v(k+1)=vt
      go to 100
      end
      

c     MISCELLANOUS SMALL SUBROUTINES

      subroutine perm(ns,ntp)
      double precision rnd
      integer ntp(ns)
      do 1 n= 1,ns
         ntp(n)=n
 1    continue        
      j=ns
 11   call rrand(rnd)
      k=int(j*rnd)+1
      jx=ntp(j)
      ntp(j)=ntp(k)
      ntp(k)=jx
      j=j-1
      if(j.gt.1) go to 11
      end

      subroutine mypack(l,icat,npack)
      integer*4 icat(32),npack
      
c     icat is a binary integer with ones for categories going left
c     and zeroes for those going right.  The sub returns npack- the integer 
c     corresponding whose binary expansion is icat.
      
      npack=0
      do 10,k=1,l
         npack=npack+icat(k)*(2**(k-1))
 10   continue
      end
      
      subroutine myunpack(l,npack,icat)
      
c npack is a long integer.  The sub. returns icat, an integer of zeroes and
c ones corresponding to the coefficients in the binary expansion of npack.
      
      integer icat(32),npack
      do j=1,32
         icat(j)=0
      end do
      n=npack
      icat(1)=mod(n,2)
      do 10 k=2,l
         n=(n-icat(k-1))/2
         icat(k)=mod(n,2)
 10   continue
      end
      
      subroutine zerv(ix,m1)
      integer ix(m1)
      do 10 n=1,m1
         ix(n)=0
 10   continue
      end
      
      subroutine zervr(rx,m1)
      double precision rx(m1)
      do 10 n=1,m1
         rx(n)=0.0d0
 10   continue
      end
      
      subroutine zerm(mx,m1,m2)
      integer mx(m1,m2)
      do 10 i=1,m1
         do 20 j=1,m2
            mx(i,j)=0
 20      continue
 10   continue
      end
      
      subroutine zermr(rx,m1,m2)
      double precision rx(m1,m2)
      do 10 i=1,m1
         do 20 j=1,m2
            rx(i,j)=0.0d0
 20      continue
 10   continue
      end

      subroutine zermd(rx,m1,m2)
      double precision rx(m1,m2)
      do 10 i=1,m1
         do 20 j=1,m2
            rx(i,j)=0.0d0
 20      continue
 10   continue
      end
      
      subroutine eqm(j,k,m,n)
      integer j(m,n),k(m,n)
      do n1=1,n
         do m1=1,m
            j(m1,n1)=k(m1,n1)
         end do
      end do
      end
      

      subroutine perm1(np,ns,tp)
      implicit double precision (a-h,o-z)
      double precision rnd, tp(ns), tx
      j=np
 11   call rrand(rnd)
      k=int(j*rnd)+1
      tx=tp(j)
      tp(j)=tp(k)
      tp(k)=tx
      j=j-1
      if(j.gt.1) go to 11
      end
      
