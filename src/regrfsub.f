C     Last change:  LB   13 Mar 2002   12:01 pm
	
c       copyright 1999 by leo Breiman
c       this is free software and can be used for any purpose. 
c       It comes with no guarantee.  
	
	
	
c	SUBROUTINE BUILDTREE
	subroutine rbuildtree(x,y,yl,mdim,nls,nsample,treemap,
     1    jdex,upper,avnode,bestcrit, nodestatus,
     2    nodepop,nodestart,nrnodes,nthsize,rsnodecost,
     3    ncase,parent,ut,v,xt,mtry,ip,
     4    mbest,cat,tgini, mind)
	
      implicit double precision (a-h,o-z)
      integer treemap(2,nrnodes),parent(nrnodes),
     $     nodestatus(nrnodes),ip(mdim),nodepop(nrnodes),
     $     nodestart(nrnodes),jdex(nsample),ncase(nsample),
     $     mbest(nrnodes),cat(mdim), mind(mdim)
	
      double precision y(nsample),bestcrit(nrnodes),x(mdim,nsample),
     1     avnode(nrnodes),xt(nsample),upper(nrnodes),
     2     v(nsample),ut(nsample),rsnodecost(nrnodes),
     3     yl(nsample),tgini(mdim)
      
      call zerv(nodestatus,nrnodes)
      call zerv(nodestart,nrnodes)
      call zerv(nodepop,nrnodes)
      call zervr(avnode,nrnodes)
	
      do n=1,nsample
         ut(n)=0
         jdex(n)=n
      end do
      
      ncur=1
      nodestart(1)=1
      nodepop(1)=nls
      nodestatus(1)=2
      
      av=0.0
      ss=0.0
      do n=1,nls
         d=y(jdex(n))
         ss=ss+(n-1)*(av-d)*(av-d)/n
         av=((n-1)*av+d)/n
      end do
      avnode(1)=av
      rsnodecost(1)=ss/nls
      
c     start main loop
      
      do 30 kbuild=1,nrnodes
         if (kbuild.gt.ncur) goto 50
         if (nodestatus(kbuild).ne.2) goto 30
         
c     initialize for next call to findbestsplit
         
         ndstart=nodestart(kbuild)
         ndend=ndstart+nodepop(kbuild)-1
         nodecnt=nodepop(kbuild)
         sumnode=nodecnt*avnode(kbuild)
         jstat=0
         
         
         call rfindbestsplit(x,xt,ut,jdex,y,mdim,nsample,
     $        ndstart,ndend,msplit,decsplit,ubest,ncase,ndendl,
     $        jstat,v,mtry,ip,sumnode,nodecnt,yl,cat, mind)

c         iprint=0
c         if(iprint.eq.1) then
c            write(*,*) kbuild,decsplit,ubest
c            write(*,*) ndstart,ndendl,ndend
c            write(*,*) ""
c         end if
         
         
         if (jstat.eq.1) then
            nodestatus(kbuild)=-1
            go to 30
         else
            mbest(kbuild)=msplit
            upper(kbuild)=ubest
            bestcrit(kbuild)=decsplit
         end if
         
         tgini(msplit)=tgini(msplit)+decsplit
                  
c     leftnode no.= ncur+1, rightnode no. = ncur+2.
         
         nodepop(ncur+1)=ndendl-ndstart+1
         nodepop(ncur+2)=ndend-ndendl
         nodestart(ncur+1)=ndstart
         nodestart(ncur+2)=ndendl+1
         
         av=0.0
         ss=0.0
         do n=ndstart,ndendl
            d=y(jdex(n))
            k=n-ndstart
            ss=ss+k*(av-d)*(av-d)/(k+1)
            av=(k*av+d)/(k+1)
         end do
         avnode(ncur+1)=av
         rsnodecost(ncur+1)=ss/nls
         
         av=0.0
         ss=0.0
         do n=ndendl+1,ndend
            d=y(jdex(n))
            k=n-ndendl-1
            ss=ss+k*(av-d)*(av-d)/(k+1)
            av=(k*av+d)/(k+1)
         end do
         avnode(ncur+2)=av
         rsnodecost(ncur+2)=ss/nls
         
c     check on nodestatus
         
         nodestatus(ncur+1)=2
         nodestatus(ncur+2)=2
         if (nodepop(ncur+1).le.nthsize)
     $        nodestatus(ncur+1)=-1
         if (nodepop(ncur+2).le.nthsize)
     $        nodestatus(ncur+2)=-1
         
         treemap(1,kbuild)=ncur+1
         treemap(2,kbuild)=ncur+2
         parent(ncur+1)=kbuild
         parent(ncur+2)=kbuild
         nodestatus(kbuild)=1
         ncur=ncur+2
         
         if (ncur.ge.nrnodes) goto 50
         
 30   continue
 50   continue

      return
      end
      
      
c	SUBROUTINE FINDBESTSPLIT
      
      subroutine rfindbestsplit(x,xt,ut,jdex,y,mdim,
     $     nsample,ndstart,ndend,msplit,decsplit,ubest,
     $     ncase,ndendl,jstat,v,mtry,ip,
     $     sumnode,nodecnt,yl,cat, mind)

      implicit double precision (a-h,o-z)
      integer ncase(nsample),jdex(nsample),ip(mdim),
     $     ncat(32),icat(32),cat(mdim), mind(mdim)
      
      double precision x(mdim,nsample),ut(nsample),xt(nsample),
     $     v(nsample),y(nsample),yl(nsample),
     $     sumcat(32),avcat(32),tavcat(32),rnd,ubestt
      
	
c       START BIG LOOP
      
      critmax=0
      ubestt = 0.0
c 200  call zerv(ip,mdim)
      non=0
      do k = 1, mdim
	 mind(k) = k
      end do

      nn = mdim
      do mt=1,mtry
         critvar=0
 100     call rrand(rnd)
	 j = int(rnd * nn) + 1
	 kv = mind(j)
	 mind(j) = mind(nn)
	 nn = nn - 1
c	 if(ip(kv).eq.1) goto 100
c	 ip(kv) = 1
	 lc=cat(kv)
         if(lc.eq.1) then
            do n=ndstart,ndend
               xt(n)=x(kv,jdex(n))
               yl(n)=y(jdex(n))
            end do
            
         else
            
            call zervr(sumcat,32)
            call zerv(ncat,32)
            do n=ndstart,ndend
               l=nint(x(kv,jdex(n)))
               d=y(jdex(n))
               sumcat(l)=sumcat(l)+d
               ncat(l)=ncat(l)+1
            end do
            do  j=1,lc
               if(ncat(j).gt.0) then
                  avcat(j)=sumcat(j)/ncat(j)
               else
                  avcat(j)=0
               end if
            end do
            do n=1,nsample
               xt(n)=avcat(nint(x(kv,jdex(n))))
               yl(n)=y(jdex(n))
            end do
         end if
         
         do n=ndstart,ndend
            v(n)=xt(n)
         end do
         do n=1,nsample
            ncase(n)=n
         end do
         call quicksort (v,ncase,ndstart,ndend,nsample)
         if(v(ndstart).ge.v(ndend)) then
            non=non+1
            if(non.ge.3*mdim) then
               jstat=1
               return
            end if
c            goto 100
	    continue
         end if
         
c     ncase(n)=case number of v nth from bottom
         
         suml=0
         sumr=sumnode
         npopl=0
         npopr=nodecnt
         
	 do nsp=ndstart,ndend-1
            d=yl(ncase(nsp))
            suml=suml+d
            sumr=sumr-d
            npopl=npopl+1
            npopr=npopr-1
            if (v(nsp).lt.v(nsp+1)) then
               crit=(suml*suml/npopl)+(sumr*sumr/npopr)
               if (crit.gt.critvar) then
                  ubestt=(v(nsp)+v(nsp+1))/2.0
                  critvar=crit
                  nbestt=nsp
               endif
            end if
         end do
         
         if(critvar.gt.critmax) then
            ubest=ubestt
            nbest=nbestt
            msplit=kv
            critmax=critvar
            do n=ndstart,ndend
               ut(n)=xt(n)
            end do
            if (cat(kv).gt.1) then
               ic=cat(kv)
               do j=1,ic
                  tavcat(j)=avcat(j)
               end do
            end if
         end if
      end do
      
      nl=ndstart-1
      do nsp=ndstart,ndend
         if(ut(nsp).le.ubest) then
            nl=nl+1
            ncase(nl)=jdex(nsp)
         end if
      end do
      
      ndendl=max0(nl,ndstart+1)
      nr=ndendl
      do nsp=ndstart,ndend
         if(ut(nsp).gt.ubest) then
            nr=nr+1
            if(nr.gt.nsample) goto 765
            ncase(nr)=jdex(nsp)
         end if
      end do
 765  continue
      
      if(ndendl.ge.ndend) ndendl=ndend-1
      
      
      do n=ndstart,ndend
         jdex(n)=ncase(n)
      end do
      lc=cat(msplit)
      if(lc.gt.1) then
         do j=1,lc
            if(tavcat(j).lt.ubest) then
               icat(j)=1
            else
               icat(j)=0
            end if
         end do
         call mypack(lc,icat,nubest)
         ubest=dble(nubest)
      end if
      
      decsplit=critmax-(sumnode*sumnode/nodecnt)
      return
      end


      subroutine rtestreebag(x,nsample,mdim,treemap,nodestatus,
     1    nrnodes,ndbigtree,ytree,upper,avnode,mbest,cat,nodex)

      implicit double precision (a-h,o-z)
      double precision x(mdim,nsample),
     1    upper(nrnodes),avnode(nrnodes),ytree(nsample)

      integer treemap(2,nrnodes),nodestatus(nrnodes),
     1    mbest(nrnodes),cat(mdim),icat(32), nodex(nsample)

      do n = 1, nsample
         nodex(n) = 0
      end do

	do n=1,nsample
	   kt=1
	   do k=1,ndbigtree
	      if(nodestatus(kt).eq.-1) then
		 ytree(n)=avnode(kt)
		 nodex(n)=kt
		 goto 100
	      end if
	      m=mbest(kt)
	      lc=cat(m)
	      if(lc.eq.1) then
		 if (x(m,n).le.upper(kt)) then 
		    kt=treemap(1,kt)
		 else
		    kt=treemap(2,kt)
		 endif 
	      else
		 mm=nint(upper(kt))
		 call myunpack(lc,mm,icat)
		 j=nint(x(m,n))
		 if(icat(j).eq.1) then 
		    kt=treemap(1,kt)
		 else
		    kt=treemap(2,kt)
		 endif
	      end if
	   end do
 100	   continue
	end do
	return
	end
	
