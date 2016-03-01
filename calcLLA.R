# Estimate derivatives using LLA (local linear approximation)
calcLLA <- function(pos,tau,delta){
	lag = shift(pos,-tau)
	lead = shift(pos,tau)
	vel = (lead - lag)/(2*tau*delta)
	acc = (lead - 2*pos + lag)/((tau^2)*(delta^2))
	posC = pos^3
	velC = vel^3
	pXvel2 = pos*(vel^2)
	vXpos2 = vel*(pos^2)
return(as.data.frame(cbind(pos,vel,acc,posC,velC,pXvel2,vXpos2)))
}