c     Copyright (C) 2001  Leo Breiman and Adele Cutler
c     This program is free software; you can redistribute it and/or
c     modify it under the terms of the GNU General Public License
c     as published by the Free Software Foundation; either version 2
c     of the License, or (at your option) any later version.

c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.
c     
c     Modified by Andy Liaw and Matt Wiener:
c     The main program is re-written as a C function to be called from R.
c     All calls to the uniform RNG is replaced with R's RNG.  Some subroutines
c     not called are excluded.  Variables and arrays declared as double as 
c     needed.  Unused variables are deleted.
c     
c     SUBROUTINE BUILDTREE
      
      subroutine buildtree(a,b,cl,cat,mdim,nsample,nclass,treemap,
     1     bestvar,bestsplit,bestsplitnext,tgini, nodestatus,nodepop,
     1     nodestart,classpop,tclasspop,tclasscat,ta,nrnodes,
     1     idmove,ndsize,ncase,parent,jin,mtry,iv,nodeclass,ndbigtree,
     1     win,wr,wc,wl,mred,nuse,mind)
      

c     Buildtree consists of repeated calls to two subroutines, Findbestsplit
c     and Movedata.  Findbestsplit does just that--it finds the best split of
c     the current node.  Movedata moves the data in the split node right and
c     left so that the data corresponding to each child node is contiguous.
c     The buildtree bookkeeping is different from that in Friedman's original
c     CART program.  ncur is the total number of nodes to date.
c     nodestatus(k)=1 if the kth node has been split.  nodestatus(k)=2 if the
c     node exists but has not yet been split, and =-1 of the node is terminal.
c     A node is terminal if its size is below a threshold value, or if it is
c     all one class, or if all the x-values are equal.  If the current node k
c     is split, then its children are numbered ncur+1 (left), and
c     ncur+2(right), ncur increases to ncur+2 and the next node to be split is
c     numbered k+1.  When no more nodes can be split, buildtree returns to the
c     main program.
      
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
     1     wl(nclass),tgini(mdim), xrand
      integer msplit
      
      msplit = 0
      call zerv(nodestatus,nrnodes)
      call zerv(nodestart,nrnodes)
      call zerv(nodepop,nrnodes)
      call zermr(classpop,nclass,nrnodes)      
      
      do 20 j=1,nclass
         classpop(j,1) = tclasspop(j)
 20   continue

      ncur = 1
      nodestart(1) = 1
      nodepop(1) = nuse
      nodestatus(1) = 2
      
c     start main loop

      do 30 kbuild = 1, nrnodes

         if (kbuild .gt. ncur) goto 50
         if (nodestatus(kbuild) .ne. 2) goto 30
         
c     initialize for next call to findbestsplit

         ndstart = nodestart(kbuild)
         ndend = ndstart + nodepop(kbuild) - 1
         do 40 j = 1, nclass
            tclasspop(j) = classpop(j,kbuild)
 40      continue
         jstat = 0

         call findbestsplit(a,b,cl,mdim,nsample,nclass,cat,ndstart,
     1        ndend,tclasspop,tclasscat,msplit,decsplit,nbest,ncase,
     1        jstat,jin,mtry,win,wr,wc,wl,mred,kbuild,mind)
         
         
         if (jstat .eq. 1) then
            nodestatus(kbuild) = -1
            goto 30
         else
            bestvar(kbuild) = msplit
            iv(msplit) = 1
            tgini(msplit) = decsplit + tgini(msplit)
            if (cat(msplit) .eq. 1) then
               bestsplit(kbuild) = a(msplit,nbest)
               bestsplitnext(kbuild) = a(msplit,nbest+1)
            else
               bestsplit(kbuild) = nbest
               bestsplitnext(kbuild) = 0
            endif
         endif
         
         call movedata(a,ta,mdim,nsample,ndstart,ndend,idmove,ncase,
     1        msplit,cat,nbest,ndendl)
         
c     leftnode no.= ncur+1, rightnode no. = ncur+2.

         nodepop(ncur+1) = ndendl - ndstart + 1
         nodepop(ncur+2) = ndend - ndendl
         nodestart(ncur+1) = ndstart
         nodestart(ncur+2) = ndendl + 1

c     find class populations in both nodes
         
         do 60 n=ndstart,ndendl
            nc = ncase(n)
            j=cl(nc)
            classpop(j,ncur+1) = classpop(j,ncur+1) + win(nc)
 60      continue
         do 70 n=ndendl+1,ndend
            nc = ncase(n)
            j = cl(nc)
            classpop(j,ncur+2) = classpop(j,ncur+2) + win(nc)
 70      continue


c     check on nodestatus

         nodestatus(ncur+1) = 2
         nodestatus(ncur+2) = 2
         if (nodepop(ncur+1).le.ndsize) nodestatus(ncur+1) = -1
         if (nodepop(ncur+2).le.ndsize) nodestatus(ncur+2) = -1
         popt1=0
         popt2=0
         do j=1,nclass
            popt1 = popt1 + classpop(j,ncur+1)
            popt2 = popt2 + classpop(j,ncur+2)
         end do
         
         do j=1,nclass
            if (classpop(j,ncur+1).eq.popt1) nodestatus(ncur+1) = -1
            if (classpop(j,ncur+2).eq.popt2) nodestatus(ncur+2) = -1
         end do

         treemap(1,kbuild) = ncur + 1
         treemap(2,kbuild) = ncur + 2
         parent(ncur+1) = kbuild
         parent(ncur+2) = kbuild
         nodestatus(kbuild) = 1
         ncur=ncur+2
         if (ncur.ge.nrnodes) goto 50
         
 30   continue
 50   continue

      ndbigtree = nrnodes
      do k=nrnodes, 1, -1
         if (nodestatus(k) .eq. 0) ndbigtree = ndbigtree - 1
         if (nodestatus(k) .eq. 2) nodestatus(k) = -1
      end do

      
      do kn = 1, ndbigtree
         if(nodestatus(kn) .eq. -1) then
            pp = 0
            do j = 1, nclass
               if (classpop(j,kn) .gt. pp) then
                  nodeclass(kn)=j
                  pp=classpop(j,kn)
               end if
c     
c     Break ties at random:
c     
               if (classpop(j,kn) .eq. pp) then
                  call rrand(xrand)
                  if (xrand .gt. 0.5) then
                     nodeclass(kn)=j
                     pp=classpop(j,kn)
                  end if
               end if

            end do
         end if
      end do
      
      end


c     SUBROUTINE FINDBESTSPLIT

c     For the best split, msplit is the variable split on. decsplit is the
c     dec. in impurity.  If msplit is numerical, nsplit is the case number
c     of value of msplit split on, and nsplitnext is the case number of the
c     next larger value of msplit.  If msplit is categorical, then nsplit is
c     the coding into an integer of the categories going left.

      subroutine findbestsplit(a,b,cl,mdim,nsample,nclass,cat,
     1     ndstart,ndend,tclasspop,tclasscat,msplit,decsplit,nbest,
     1     ncase,jstat,jin,mtry,win,wr,wc,wl,mred,kbuild,mind)
      implicit double precision(a-h,o-z)      
      integer a(mdim,nsample),cl(nsample),cat(mdim),
     1     ncase(nsample),b(mdim,nsample),jin(nsample), nn, j          
      double precision tclasspop(nclass),tclasscat(nclass,32),
     1     win(nsample),
     1     wr(nclass),wc(nclass),wl(nclass), xrand
      integer mind(mred)
c     real pno, pdo, rrd, rld
      
c     compute initial values of numerator and denominator of Gini
      
      pno = 0.0
      pdo = 0.0
      do 10 j = 1, nclass
         pno = pno + tclasspop(j) * tclasspop(j)
         pdo = pdo + tclasspop(j)
 10   continue
      crit0 = pno / pdo
      jstat = 0
c     zz=rrand()
c     call rrand(zz)
      
c     start main loop through variables to find best split
      
      critmax = -1.0e20
      
      do k = 1, mred
         mind(k) = k
      end do

      nn = mred
      do 20 mt = 1, mtry
 200     continue
c     
c     sampling mtry variables w/o replacement.
c     
         call rrand(xrand)
         j = int(nn * xrand) + 1
         mvar = mind(j)
         mind(j) = mind(nn)
         mind(nn) = mvar
         nn = nn - 1
c     
         if(cat(mvar) .eq. 1) then
            rrn = pno
            rrd = pdo
            rln = 0
            rld = 0
            call zervr(wl, nclass)
            do 50 j = 1, nclass
               wr(j) = tclasspop(j)
 50         continue
            critvar = -1e20
            
            do 60 nsp = ndstart, ndend-1
               nc=a(mvar, nsp)
               u = win(nc)
               k = cl(nc)
               rln = rln + u * (2 * wl(k) + u)
               rrn = rrn + u * (-2 * wr(k) + u)
               rld = rld + u
               rrd = rrd - u
               wl(k) = wl(k) + u
               wr(k) = wr(k) - u
               
               if (b(mvar, nc) .lt. b(mvar,a(mvar, nsp + 1))) then
                  if(dmin1(rrd, rld) .gt. 1.0e-5) then
                     crit = (rln / rld) + (rrn / rrd)
                     if (crit .gt. critvar) then
                        nbestvar=nsp
                        critvar=crit
                     endif
c     
c     Break ties at random:
c     
                     if (crit .eq. critvar) then
                        call rrand(xrand)
                        if (xrand .gt. 0.5) then
                           nbestvar = nsp
                           critvar = crit
                        end if
                     end if
c     
                  end if
               end if
 60         continue
 65         continue

            if (critvar .gt. critmax) then
               msplit = mvar
               nbest = nbestvar
               critmax = critvar
            end if
c     
c     Break ties at random:
c     
            if (critvar .eq. critmax) then
               call rrand(xrand)
               if (xrand .gt. 0.5) then
                  msplit = mvar
                  nbest = nbestvar
                  critmax = critvar
               end if
            end if

         else

c     compute the decrease in impurity given by categorical splits

            lcat = cat(mvar)
            call zermr(tclasscat, nclass,32)
            do 70 nsp = ndstart, ndend
               nc = ncase(nsp)
               l = a(mvar,ncase(nsp))
               tclasscat(cl(nc), l) = tclasscat(cl(nc), l) + win(nc)
 70         continue
            nnz = 0
            do i = 1, lcat
               su = 0
               do j = 1, nclass
                  su = su + tclasscat(j, i)
               end do
               if(su .gt. 0) nnz = nnz + 1
            end do
            if (nnz.eq.1) then
               critvar = -1.0e25
            else
               call catmax(pno, pdo, tclasscat, tclasspop, nclass, lcat,
     1              nbestvar, critvar)
            end if
            
c     this last subroutine returns those categories going left in the best split. 
c     This is coded into a long integer (see under subroutine catmax below for 
c     details). 
            
            if (critvar .gt. critmax) then
               msplit = mvar
               nbest = nbestvar
               critmax = critvar
            endif
c     
c     Break ties at random:
c     
            if (critvar .eq. critmax) then
               call rrand(xrand)
               if (xrand .gt. 0.5) then
                  msplit = mvar
                  nbest = nbestvar
                  critmax = critvar
               end if
            end if
         end if
 20   continue
 25   continue                
      decsplit = critmax - crit0
      if (critmax .lt. -1.0e10) jstat = 1

      return
      end
      
C     SUBROUTINE CATMAX

      subroutine catmax(pno, pdo, tclasscat, tclasspop, nclass, lcat,
     1     ncatsplit, rmaxdec)
      
c     this subroutine finds the best categorical split of a categorical variable
c     with lcat categories, nclass classes and tclasscat(j,l) is the number of 
c     cases in class j with category value l. The method used is an exhaustive 
c     search over all partitions of the category values.  For the two class 
c     problem, there is a faster exact algorithm we will add later.  
c     If lcat.ge.10, the exhaustive search gets slow and there is a faster 
c     iterative algorithm we can add later.
      
      parameter(jmax=1000)
      implicit double precision(a-h,o-z)
      double precision tclasscat(nclass,32),tclasspop(nclass),
     1     tmpclass(jmax), xrand
      integer icat(32) 
      
      rmaxdec = -1e20
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
         if (tdec .gt. rmaxdec) then
            rmaxdec=tdec
            ncatsplit=n
         endif
c     
c     Break ties at random:
c     
         if (tdec .eq. rmaxdec) then
            call rrand(xrand)
            if (xrand .gt. 0.5) then
               rmaxdec=tdec
               ncatsplit=n
            end if
         end if

 10   continue
      
      return
      end

c     
c     -------------------------------------------------------
      subroutine catmaxb(tclasscat,tclasspop,xc,cp,cm,kcat,
     1     nclass,lcat,maxcat,ncatsplit,critmax,pdo,nhit,dn)
      implicit double precision (a-h,o-z)
      double precision tclasscat(nclass,maxcat),tclasspop(nclass),
     $     xc(maxcat),cp(maxcat),cm(maxcat),dn(maxcat)
      double precision critmax, pdo
      integer kcat(maxcat),ncatsplit(maxcat)
      integer nclass, lcat, maxcat, nhit, l, nk, j, n, i, k
      double precision bestsplit, rrd, rld, rln, rrn, crit
      
      nhit = 0
      do l = 1, lcat
         if (dn(l) .gt. 0) then
            xc(l) = tclasscat(1,l) / dn(l)
         else
            xc(l)=0
         end if
      end do
      do nk = 1, lcat
         kcat(nk) = nk
      enddo
      call quicksort(xc, kcat, 1, lcat, maxcat)
      do j = 1, nclass
         cp(j) = 0
         cm(j) = tclasspop(j)
      end do
      rrd = pdo
      rld = 0
      do n = 1, lcat - 1
         rld = rld + dn(kcat(n))
         rrd = rrd - dn(kcat(n))
         do k = 1, nclass
            cp(k) = cp(k) + tclasscat(k,kcat(n))
            cm(k) = cm(k) - tclasscat(k,kcat(n))
         end do
         rln = 0
         rrn = 0
         do k = 1, nclass
            rln = rln + cp(k)**2
            rrn = rrn + cm(k)**2
         end do
         if (xc(n) .lt. xc(n+1)) then
            if(dmin1(rrd, rld) .gt. 1) then
               crit = (rln / rld) + (rrn / rrd)
               if (crit .gt. critmax) then
                  critmax = crit
                  bestsplit = .5*(xc(n) + xc(n+1))
                  nhit = 1
               end if
            end if
         end if
      end do			!n
      if (nhit.eq.1) then
         call zerv(ncatsplit,maxcat)
         do l = 1, lcat
            if (dn(l) .gt. 0) then
               xc(l) = tclasscat(1,l) / dn(l)
            else
               xc(l) = 0
            end if
         end do
         do i = 1, lcat
            if(xc(i) .lt. bestsplit) ncatsplit(i) = 1
         end do
      end if
      return
      end
      
c     -------------------------------------------------------
      subroutine catmaxap(tclasscat,kcat,
     $     nclass,lcat,maxcat,critmax,dn,nhit)
c     
c     this returns those categories going left in the best split. 
c     This is coded into a vector of length maxcat (see under catmax 
c     below for details). 
c     
      implicit double precision (a-h, o-z)
      double precision tclasscat(nclass,maxcat),dn(maxcat),ncr0(100),
     $     ncl(100),ncr(100),ncl0(100),ntr,ntl,ntl0,ntr0
c     
      integer kcat(maxcat)
c     
      double precision critmax, crit, critnew
      integer nclass, lcat, maxcat, nhit, l, j, k, nchange
      double precision xrand
      
      nhit = 0
      do l = 1, lcat
         call rrand(xrand)
         if (xrand .le. 0.5) then
            kcat(l) = -1
         else
            kcat(l) = 1
         end if
      end do
      do j = 1, nclass
         ncl(j) = 0
         ncr(j) = 0
         ntl = 0
         ntr = 0
      end do
      do l = 1, lcat
         do j = 1, nclass
            if (kcat(l) .eq. -1) then
               ncl(j) = ncl(j) + tclasscat(j,l)
               ntl = ntl + dn(l)
            else
               ncr(j) = ncr(j) + tclasscat(j,l)
               ntr = ntr + dn(l)
            end if
         end do
      end do
      crit = 0
      do j = 1, nclass
         crit = crit + (ncl(j) * ncl(j) / ntl) + (ncr(j) * ncr(j) / ntr)
      end do
      do k = 1, 1000
         nchange = 0
         do l = 1, lcat
            do j = 1, nclass
               ncl0(j) = ncl(j)
               ncr0(j) = ncr(j)
               ntl0 = ntl
               ntr0 = ntr
               ncl(j) = ncl(j) + kcat(l) * tclasscat(j,l)
               ncr(j) = ncr(j) - kcat(l) * tclasscat(j,l)
               ntl = ntl + kcat(l) * dn(l)
               ntr = ntr - kcat(l) * dn(l)
            end do
            critnew = 0
            do j = 1, nclass
               critnew = critnew + (ncl(j) * ncl(j) / ntl) + 
     1              (ncr(j) * ncr(j) / ntr)
            end do
            if(critnew .gt. critmax) then
               critmax = critnew
               nchange = nchange + 1
               kcat(l) = -kcat(l)
               nhit = 1
            else
               do j = 1, nclass
                  ncl(j) = ncl0(j)
                  ncr(j) = ncr0(j)
                  ntl = ntl0
                  ntr = ntr0
               end do
            end if
         end do                 !l
         if (nchange .eq. 0) goto 101
      end do			!k
 101  continue
      return
      end
      
c     SUBROUTINE MOVEDATA     
      
c     This subroutine is the heart of the buildtree construction. Based on the 
c     best split the data in the part of the a matrix corresponding to the 
c     current node is moved to the left if it belongs to the left child and 
c     right if it belongs to the right child.
      
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
      
c     npack is a long integer.  The sub. returns icat, an integer of zeroes and
c     ones corresponding to the coefficients in the binary expansion of npack.
      
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
