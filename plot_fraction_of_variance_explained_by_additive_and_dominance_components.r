plot_variance_contributions <- function(freq, a, d, lwd=4)
{
	p <- freq
	q <- 1 - freq 
	V_A <- 2 * p * q * ( a + (q - p) * d)^2
	V_D <- (2 * p * q * d)^2
	V <- (V_A + V_D)

	print(V_A/V)
	print(V_D/V)

	par(mfrow=c(1,2))
	# Proportions
	plot(p, V_A/V, type='l', col='cornflowerblue', ylim=c(0,1),
		ylab = 'fraction of total variance',
		xlab='allele frequency', lwd=lwd)
	points(p, V_D/V, type='l', col='indianred3', lwd=lwd)
	legend(x='topleft',
		col=c('cornflowerblue', 'indianred3'),
		legend=c('$\\frac{V_{A}}{V}$', '$\\frac{V_{D}}{V}$'),
           lty=c(1,1,2,1,2), lwd=lwd, bty='n', y.intersp=1.5)

	# Contributions
	plot(p, V, type='l',
		ylab = 'fraction of variance explained',
		xlab='allele frequency', lwd=lwd)
	points(p, V_A, col='cornflowerblue', type='l', lwd=lwd)
	points(p, V_D, type='l', col='indianred3', lwd=lwd)

}

freq <- seq(0,1, length.out = 100)
plot_variance_contributions(freq, 0, 1)

freq <- seq(0,1, length.out = 100)
plot_variance_contributions(freq, 1, 0.5)

freq <- seq(0,1, length.out = 100)
plot_variance_contributions(freq, 1, -0.5)

freq <- c(0.3,0.5)
plot_variance_contributions(freq, 0, 1)

freq <- c(0.3, 0.5)
plot_variance_contributions(freq, 1, 0.5)
