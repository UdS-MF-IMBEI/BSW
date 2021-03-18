df <- data.frame(y = rep(c(0, 1), each = 250), 
                 x = rep(c(0, 1, 0, 1), times = c(200, 50, 50, 200))
)
RR <- (200 * 250) / (50 * 250)
SE <- sqrt((1/200 + 1/50) - (1/250 + 1/250))
fit <- bsw(y ~ x, df)
out <- summary(fit)

test_that(desc = "Estimated relative risk is equal to 4",
          code = {
            expect_equal(object = unname(exp(coef(fit)[2])),
                         expected = RR)
            expect_equal(object = unname(out$std.err[2]), 
                         expected = SE)
            expect_equal(object = unname(out$z.value[2]), 
                         expected = log(RR) / SE)
            expect_equal(object = unname(exp(confint(fit)[2,])), 
                         expected = exp(log(RR) + SE * qnorm(c(0.025, 0.975))))
          }
)