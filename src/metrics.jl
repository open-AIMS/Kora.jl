RMSE(y_hat, y) = sqrt(mean((y_hat .- y).^2))
R2(y_hat, y) = 1 - (sum((y .- y_hat).^2) / sum((y .- mean(y)).^2))