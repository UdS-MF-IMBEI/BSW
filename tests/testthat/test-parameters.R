# Manual computation
df <- data.frame(y = rep(c(0, 1), each = 250),
                 x = rep(c(0, 1, 0, 1), times = c(200, 50, 50, 200)))
RR <- (200 * 250) / (50 * 250)
names(RR) <- "x"
SE <- sqrt((1/200 + 1/50) - (1/250 + 1/250))
names(SE) <- "x"
Z <- log(RR) / SE
CI <- exp(log(RR) + SE * qnorm(c(0.025, 0.975)))
names(CI) <- c("2.5%", "97.5%")

# Algorithmic computation
fit <- bsw(y ~ x, df)
cf <- coef(fit)
rr <- exp(cf[2])
se <- sqrt(diag(solve(hess(cf, fit@y, fit@x))))[2]
z <- cf[2] / se
ci <- exp(confint(fit)[2,])

# Check numerical equality
test_that(desc = "Check numerical equality",
          code = {
            expect_equal(object = rr,
                         expected = RR)
            expect_equal(object = se,
                         expected = SE)
            expect_equal(object = z,
                         expected = Z)
            expect_equal(object = ci,
                         expected = CI)
          }
)
