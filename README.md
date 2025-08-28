# Denmark-Fertility-Forecasting-

This project forecasts monthly births in Denmark for 2025 using official data from 2010 to 2025 with Seasonal ARIMA (SARIMA) and dynamic models with seasonal regressors (SARIMAX), incorporating external shocks and structural breaks linked to the 2014 “Do it for Denmark” fertility campaign and the post-COVID-19 effect. SARIMAX with campaign and COVID dummies achieved the best performance, showing that accounting for structural breaks improves demographic forecasts. The SARIMAX(1,1,1)(0,1,1)[12] model, including both shock regressors, produced the lowest error metrics and the best BIC and AIC values.

📄 [Read the full paper (PDF)](./report.pdf)  
📘 [Open the R Studio file (R)](./code.R)
