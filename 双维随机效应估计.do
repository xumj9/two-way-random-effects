
//说明：以陈强书traffic.dta的数据为例，解释变量为X中的变量，被解释变量为Y中的变量，betahat为估计系数矩阵，varbetahat为估计标准差
sort state year
putmata X=(unrate spircons beertax)
putmata Y=fatal
xtset state year
scalar T=r(tmax)-r(tmin)
egen tmin=min(year)
egen tmax=max(year)
gen t=tmax-tmin+1
putmata T=t

mata
	 obs=rows(Y)
	 T=T[1]
	 n=obs/T
	 K=cols(X)
	 in=J(1,n,1)
	 it=J(1,T,1)
	 
	 Q1=(I(n)-cross(in,in)/n)#(I(T)-cross(it,it)/T)
	 Q2=(I(n)-cross(in,in)/n)#(cross(it,it)/T)
	 Q3=(cross(in,in)/n)#(I(T)-cross(it,it)/T)
	 Q4=(cross(in,in)/n)#(cross(it,it)/T)
	 
	 defla1=(n-1)*(T-1)-K
	 Lamda1=(Y'*Q1*Y-Y'*Q1*X*cholinv(X'*Q1*X)*X'*Q1*Y)/defla1 
	 defla2=(n-1)-K
	 Lamda2=(Y'*Q2*Y-Y'*Q2*X*cholinv(X'*Q2*X)*X'*Q2*Y)/defla2 
	 defla3=(T-1)-K
	 Lamda3=(Y'*Q3*Y-Y'*Q3*X*cholinv(X'*Q3*X)*X'*Q3*Y)/defla3 
	 Lamda4=Lamda3+Lamda2-Lamda1
	 
	 Vminus=Q1/Lamda1+Q2/Lamda2+Q3/Lamda3+Q4/Lamda4

	 betahat=cholinv(X'*Vminus*X)*X'*Vminus*Y
	 betahat
	 varbetahat=diag(Lamda1*cholinv(X'*Vminus*X))
	 varbetahat
end
