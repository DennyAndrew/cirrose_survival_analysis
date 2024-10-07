# Base de dados extraida do kaglle referente 
# Projeto de Pesquisa

# Não esqueça de conferir seu diretório atual.
# libs necessárias
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(gamlss)
library(gamlss.cens)

# Devido a algun processos dos pacotes, precisamos definir uma semente, para 
# assegurar reprodutividade
set.seed(007)

# Base de dados extraida do kaglle referente 
base <- read.csv("cirrhosis.csv", h = T)

# Um total de 424 pacientes com PBC, encaminhados à Clínica Mayo durante esse 
# intervalo de dez anos, preencheram os critérios de elegibilidade para o ensaio
# randomizado controlado por placebo do medicamento D-penicilamina. Os primeiros
# 312 casos do conjunto de dados participaram do ensaio randomizado e contêm 
# dados praticamente completos. Os 112 casos adicionais não participaram do 
# ensaio clínico, mas consentiram que as medidas básicas fossem registradas e 
# acompanhadas para sobrevivência. Seis desses casos foram perdidos no 
# acompanhamento logo após o diagnóstico, portanto, os dados aqui referem-se a 
# mais 106 casos, bem como aos 312 participantes randomizados.

# Variáveis:
## ID: identificador
## N_days: número de dias entre o registro e o que ocorrer primeiro entre morte, 
# transplante ou análise do estudo em julho de 1986
## Status: status do paciente C (censurado), CL (censurado devido ao tx do fígado) ou D (óbito)
## Drug: tipo de medicamento D-penicilamina ou placebo
## Age: idade em dias
## Sex: sexo do paciente
## Ascites: presença de ascite N (Não) ou S (Sim)
## Hepatomegaly: presença de hepatomegalia N (Não) ou Y (Sim)
## Spiders: presença de aranhas N (Não) ou Y (Sim)
## Edema: presença de edema N (sem edema e sem terapia diurética para edema),
# S (edema presente sem diuréticos ou edema resolvido com diuréticos) ou Y (edema apesar da terapia diurética)
## Bilirubin: bilirrubina sérica em [mg/dl]
## Cholesterol: colesterol sérico em [mg/dl]
## Albumin: albumina em [gm/dl]
## Copper: cobre na urina em [ug/dia]
## Alk_Phos: fosfatase alcalina em [U/litro]
## SGOT: Sem definição
## Tryglicerides: triglicerídeos em [mg/dl]
## Platelets: plaquetas por cúbico [ml/1000]
## Prothrombin: tempo de protrombina em segundos [s]
## Stage: estágio histológico da doença (1, 2, 3 ou 4)

# Vamos inicialmente retirar a variável de ID, pois não fornece informações
base <- base[, !names(base) %in% "ID"]

# Para a variável Drug os casos NA não nos trazem nenhum tipo de informações e
# referen-se aos casos que não participaram do estudo clinico, logo essas observações
# serão retiradas da nossa base, mantendo somente os clientes que receberam a droga
# ou placebo
base <- base %>%
  filter(!is.na(base[['Drug']]))
# O evento de interesse trata-se do óbito do paciente, portanto os casos em 
# que ocorreram censura por transplante, vamos tratar todos como censura, pois 
# o tempo de evento ainda seria superior ao tempo de censura.
base <- base %>%
  mutate(Status = ifelse(Status %in% c("CL", "C"), "C", "E"))
base$Cens <- ifelse(base$Status == "E", 1, 0)
base <- base[, !names(base) %in% "Status"]

# padronizando a variável idade para anos
base$Age <- round(base$Age/365, 0)

median(base$Age, na.rm = TRUE)
# temos que a mediana é de 50, sendo assim vamos categorizar em dois grupos;
# < 50 anos e > 50 anos
base$Age <- ifelse(base$Age < 50, "<50anos", ">50anos")

# Para a variável Endema, vamos separar apenas entre presença ou não presença
base <- base %>%
  mutate(Edema = ifelse(Edema %in% c("Y", "S"), "YES", "NO"))

# Bilirubin - categorização pela mediana
median(base$Bilirubin, na.rm = TRUE)
base$Bilirubin <- ifelse(base$Bilirubin <= 1.35, "<=1.35[mg/dl]", ">1.35[mg/dl]")

# Cholesterol - categorização pela mediana (28 NAs)
median(base$Cholesterol, na.rm = TRUE)
base$Cholesterol <- ifelse(base$Cholesterol <= 309.5, "<=309.5[mg/dl]", ">309.5[mg/dl]")

# Albumin - categorização pela mediana
median(base$Albumin, na.rm = TRUE)
base$Albumin <- ifelse(base$Albumin <= 3.55, "<=3.55[mg/dl]", ">3.55[mg/dl]")

# Copper - categorização pela mediana
median(base$Copper, na.rm = TRUE)
base$Copper <- ifelse(base$Copper <= 73, "<=73[ug/day]", ">73[ug/day]")

# Alk_Phos - categorização pela mediana
median(base$Alk_Phos, na.rm = TRUE)
base$Alk_Phos <- ifelse(base$Alk_Phos <= 1259, "<=1259[U/liter]", ">1259[U/liter]")

# SGOT - categorização pela mediana
median(base$SGOT, na.rm = TRUE)
base$SGOT <- ifelse(base$SGOT <= 114.7, "<=114.7[U/ml]", ">114.7[U/ml]")

# Tryglicerides - categorização pela mediana (30 NAs)
median(base$Tryglicerides, na.rm = TRUE)
base$Tryglicerides <- ifelse(base$Tryglicerides <= 108, "<=108[mg/dl]", ">108[mg/dl]")

# Platelets - categorização pela mediana
median(base$Platelets, na.rm = TRUE)
base$Platelets <- ifelse(base$Platelets <= 257, "<=257[ml/1000]", ">257[ml/1000]")

# Prothrombin - categorização pela mediana
median(base$Prothrombin, na.rm = TRUE)
base$Prothrombin <- ifelse(base$Prothrombin <= 10.6, "<=10.6s", ">10.6s")

# Análise exploratória ------------------------------------------------------------------------------

# Freq. relativa do sexo dos pacientes
# Calcular a frequência relativa
dados_frequencia <- base %>%
  count(Sex) %>%
  mutate(frequencia_relativa = n / sum(n))

# Criando o gráfico de pizza com frequência relativa
ggplot(dados_frequencia, aes(x = "", y = frequencia_relativa, fill = Sex)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  geom_text(aes(label = scales::percent(frequencia_relativa, accuracy = 0.1)),
            position = position_stack(vjust = 0.5)) +
  labs(title = "Frequência Relativa por Sexo") +
  scale_fill_manual(values = c("M" = "#DCDCDC", "F" = "#A9A9A9")) +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))  # Centraliza o título

# Calcular a frequência relativa
# Tipo de tratamento por sexo do paciente ----------------------------------------------------------------------
dados_frequencia <- base %>%
  count(Drug, Sex) %>%
  group_by(Drug) %>%
  mutate(frequencia_relativa = n / sum(n)) %>%
  ungroup()

dados_frequencia <- dados_frequencia %>%
  mutate(frequencia_relativa_formatada = scales::percent(frequencia_relativa, accuracy = 0.1))

# Criando o gráfico de barras com frequência relativa e cores personalizadas
ggplot(dados_frequencia, aes(x = Drug, y = frequencia_relativa, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = frequencia_relativa_formatada),
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Tipo de tratamento x sexo do paciente",
       x = "Tratamento",
       y = "Frequência Relativa",
       fill = "Sex") +
  scale_fill_manual(values = c("M" = "#DCDCDC", "F" = "#A9A9A9")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5))  # Centraliza o título

# vemos que a quantidade de mulheres no estudo é significativamente superior a quantidade de 
# homens no estudo, isso pode se dar pelo fato da doença ser realmente mais comuns em mulheres.

# Tipo de tratamento por idade ----------------------------------------------------------------------
dados_frequencia <- base %>%
  count(Drug, Age) %>%
  group_by(Drug) %>%
  mutate(frequencia_relativa = n / sum(n)) %>%
  ungroup()

dados_frequencia <- dados_frequencia %>%
  mutate(frequencia_relativa_formatada = scales::percent(frequencia_relativa, accuracy = 0.1))

# Criando o gráfico de barras com frequência relativa e cores personalizadas
ggplot(dados_frequencia, aes(x = Drug, y = frequencia_relativa, fill = Age)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = frequencia_relativa_formatada),
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Tipo de tratamento x idade do paciente",
       x = "Tratamento",
       y = "Frequência Relativa",
       fill = "Idade do paciente") +
  scale_fill_manual(values = c("<50anos" = "#DCDCDC", ">50anos" = "#A9A9A9")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5))  # Centraliza o título

# o estudo foi bem balanceado em relação a distribuição dos tratamentos por idade do paciente
# Copper por sexo ----------------------------------------------------------------------
# retirando casos missing
base_copper <- base %>%
  filter(!is.na(base[['Copper']]))

dados_frequencia <- base_copper %>%
  count(Copper, Sex) %>%
  group_by(Copper) %>%
  mutate(frequencia_relativa = n / sum(n)) %>%
  ungroup()

dados_frequencia <- dados_frequencia %>%
  mutate(frequencia_relativa_formatada = scales::percent(frequencia_relativa, accuracy = 0.1))

# Criando o gráfico de barras com frequência relativa e cores personalizadas
ggplot(dados_frequencia, aes(x = Copper, y = frequencia_relativa, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = frequencia_relativa_formatada),
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Concentração de cobre na urina x sexo do paciente",
       x = "Concentração de cobre na urina",
       y = "Frequência Relativa",
       fill = "Sex") +
  scale_fill_manual(values = c("F" = "#DCDCDC", "M" = "#A9A9A9")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5))  # Centraliza o título

# fica um pouco complicado de comparar, pois a quantidade de pacientes do sexo masculino
# é muito inferior ao sexo feminino.

# Copper por tratamento ----------------------------------------------------------------------
# retirando casos missing
base_copper <- base %>%
  filter(!is.na(base[['Copper']]))

dados_frequencia <- base_copper %>%
  count(Copper, Drug) %>%
  group_by(Copper) %>%
  mutate(frequencia_relativa = n / sum(n)) %>%
  ungroup()

dados_frequencia <- dados_frequencia %>%
  mutate(frequencia_relativa_formatada = scales::percent(frequencia_relativa, accuracy = 0.1))

# Criando o gráfico de barras com frequência relativa e cores personalizadas
ggplot(dados_frequencia, aes(x = Copper, y = frequencia_relativa, fill = Drug)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = frequencia_relativa_formatada),
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Concentração de cobre na urina x Tratamento",
       x = "Concentração de cobre na urina",
       y = "Frequência Relativa",
       fill = "Tratamento") +
  scale_fill_manual(values = c("D-penicillamine" = "#DCDCDC", "Placebo" = "#A9A9A9")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5))  # Centraliza o título

# fica um pouco complicado de comparar, pois a quantidade de pacientes do sexo masculino
# é muito inferior ao sexo feminino.

## Seleção de variáveis
# para o método de seleção de variáveis, vamos ajustar estimadores não paramétricos
# de Kaplan-Meier e aplicar um teste de log rank, assim saberemos se as curvas
# de sobrevivência são significativamente diferentes, considerando os grupos de 
# cada variável.
tempos <- base$N_Days
cens <- base$Cens
splots1 <- list()
splots2 <- list()
splots3 <- list()
splots4 <- list()

# Drug (p = 0,75)
ekm1 <- survfit(Surv(tempos, cens) ~ Drug, data = base)
splots1[[1]] <- ggsurvplot(ekm1, pval = TRUE, title = "Drug")

# Age (p < 0,0001)
ekm2 <- survfit(Surv(tempos, cens) ~ Age, data = base)
splots1[[2]] <- ggsurvplot(ekm2, pval = TRUE, title = "Age")

# Sex (p = 0,039)
ekm3 <- survfit(Surv(tempos, cens) ~ Sex, data = base)
splots1[[3]] <- ggsurvplot(ekm3, pval = TRUE, title = "Sex")

# Ascites (p < 0,0001)
ekm4 <- survfit(Surv(tempos, cens) ~ Ascites, data = base)
splots1[[4]] <- ggsurvplot(ekm4, pval = TRUE, title = "Ascites")

# Hepatomegaly (p < 0,0001) 
ekm5 <- survfit(Surv(tempos, cens) ~ Hepatomegaly, data = base)
splots2[[1]] <- ggsurvplot(ekm5, pval = TRUE, title = "Hepatomegaly")

# Spiders (p < 0,0001)
ekm6 <- survfit(Surv(tempos, cens) ~ Spiders, data = base)
splots2[[2]] <- ggsurvplot(ekm6, pval = TRUE, title = "Spiders")

# Edema (p < 0,0001)
ekm7 <- survfit(Surv(tempos, cens) ~ Edema, data = base)
splots2[[3]] <- ggsurvplot(ekm7, pval = TRUE, title = "Edema")

# Bilirubin (p < 0,0001)
ekm8 <- survfit(Surv(tempos, cens) ~ Bilirubin, data = base)
splots2[[4]] <- ggsurvplot(ekm8, pval = TRUE, title = "Bilirubin")

# Cholesterol (p = 0,0077)
ekm9 <- survfit(Surv(tempos, cens) ~ Cholesterol, data = base)
splots3[[1]] <- ggsurvplot(ekm9, pval = TRUE, title = "Cholesterol")

# Albumin (p < 0,0001)
ekm10 <- survfit(Surv(tempos, cens) ~ Albumin, data = base)
splots3[[2]] <- ggsurvplot(ekm10, pval = TRUE, title = "Albumin")

# Copper (p < 0,0001)
ekm11 <- survfit(Surv(tempos, cens) ~ Copper, data = base)
splots3[[3]] <- ggsurvplot(ekm11, pval = TRUE, title = "Copper")

# Alk_Phos (p = 0,0039)
ekm12 <- survfit(Surv(tempos, cens) ~ Alk_Phos, data = base)
splots3[[4]] <- ggsurvplot(ekm12, pval = TRUE, title = "Alk_Phos")

# SGOT (p < 0,0001)
ekm13 <- survfit(Surv(tempos, cens) ~ SGOT, data = base)
splots4[[1]] <- ggsurvplot(ekm13, pval = TRUE, title = "SGOT")

# Tryglicerides (p = 0,002)
ekm14 <- survfit(Surv(tempos, cens) ~ Tryglicerides, data = base)
splots4[[2]] <- ggsurvplot(ekm14, pval = TRUE, title = "Tryglicerides")

# Platelets (p = 0,0028)
ekm15 <- survfit(Surv(tempos, cens) ~ Platelets, data = base)
splots4[[3]] <- ggsurvplot(ekm15, pval = TRUE, title = "Platelets")

# Prothrombin (p < 0,0001) 
ekm16 <- survfit(Surv(tempos, cens) ~ Prothrombin, data = base)
splots4[[4]] <- ggsurvplot(ekm16, pval = TRUE, title = "Prothrombin")

# Visualização para artigo
arrange_ggsurvplots(splots1, print = TRUE, ncol = 2, nrow = 2)
arrange_ggsurvplots(splots2, print = TRUE, ncol = 2, nrow = 2)
arrange_ggsurvplots(splots3, print = TRUE, ncol = 2, nrow = 2)
arrange_ggsurvplots(splots4, print = TRUE, ncol = 2, nrow = 2)

# De acordo com o teste de log rank, temos que apenas a variável referente a 
# droga utilizada no ensaio não teve uma diferença significativa entre as curvas
# de sobrevivência

dados <- base[, (names(base) %in% c('N_Days', 'Cens', 'Age', 'Sex', 'Ascites',
                                    'Hepatomegaly', 'Spiders', 'Edema', 'Bilirubin',
                                    'Albumin', 'Cholesterol', 'Copper', 'SGOT',
                                    'Prothrombin', 'Tryglicerides', 'Platelets'))]

# retirando as observações NA das var
dados <- dados %>%
  filter(!is.na(dados[['Copper']]))
dados <- dados %>%
  filter(!is.na(dados[['Cholesterol']]))
dados <- dados %>%
  filter(!is.na(dados[['Tryglicerides']]))
dados <- dados %>%
  filter(!is.na(dados[['Platelets']]))

# Vamos testar alguns ajustes de modelos e selecionar o que apresentar o melhor
# ajuste:
tempos <- dados$N_Days
cens <- dados$Cens


# Exponencial
Ajuste_Exp <- gamlss(Surv(tempos, cens) ~ Age+Sex+Ascites+Hepatomegaly+Spiders+Edema+Bilirubin+Albumin+Cholesterol+Copper+SGOT+Prothrombin+Tryglicerides+Platelets,
                     data = dados, family = cens(EXP),
                     method = mixed(10, 50))

# Weibull
Ajuste_Weibull <- update(Ajuste_Exp, family = cens(WEI))

# Lognormal
Ajuste_Lognormal <- update(Ajuste_Exp, family = cens(LOGNO))

# Gama
Ajuste_Gama <- update(Ajuste_Exp, family = cens(GG))

# Para selecionar o melhor ajsute, vamos realizar uma análise de resíduos
plot(Ajuste_Exp) 
plot(Ajuste_Weibull) 
plot(Ajuste_Lognormal) # Apresenta um certo nível de aderencia
plot(Ajuste_Gama) 

# podemos utilizar uma outra abordagem para análise de resíduos, com esse método
# podemos observar as bandas de confiança, permitindo uma análise mais ampla.
par(mfrow = c(1, 1))
wp(Ajuste_Exp)
wp(Ajuste_Weibull)
wp(Ajuste_Lognormal)
wp(Ajuste_Gama)

# Os ajsutes com distribuição Weibull, Log-normal e Gama apresentam um certo nível
# de qualidade. Dado isso vamos presseguir com uma análise de AIC:
aic_exp <- extractAIC(Ajuste_Exp)[2]
aic_weibull <- extractAIC(Ajuste_Weibull)[2]
aic_lognorm <- extractAIC(Ajuste_Lognormal)[2]
aic_gama <- extractAIC(Ajuste_Gama)[2]
round(cbind(aic_exp, aic_weibull, aic_lognorm, aic_gama),1)

# Os ajustes log-normal e gama apresentam um valor de aic muito semelhantes,
# sendo assim, optaremos pelo ajuste log normal, por ser o modelo mais simples, 

summary(Ajuste_Lognormal)
# Pelo summary do modelo, algumas variáveis não apresentaram significância
# perante o efeito das demais covariáveis, dessa forma será ajustado um modelo
# desconsiderando essas variáveis e comparado ao ajsute atual.

dados2 <- dados[, !names(dados) %in% "Sex"]
dados2 <- dados[, !names(dados) %in% "Hepatomegaly"]
dados2 <- dados[, !names(dados) %in% "Cholesterol"]
dados2 <- dados[, !names(dados) %in% "Tryglicerides"]
dados2 <- dados[, !names(dados) %in% "Platelets"]

# Ajuste de modelo log normal
tempos <- dados2$N_Days
cens <- dados2$Cens
Ajuste_lognormal2 <- gamlss(Surv(tempos, cens) ~ Age+Spiders+Ascites+Edema+Bilirubin+Albumin+Copper+SGOT+Prothrombin,
                            data = dados2, family = cens(LOGNO),
                            method = mixed(10, 50))
summary(Ajuste_lognormal2)

# Todas as corariáveis do ajuste apresentam significância a 5%.

# Comparando aic desse ajuste com o ajuste considerando todas as cováriaveis
aic_lognorm2 <- extractAIC(Ajuste_lognormal2)[2]
round(cbind(aic_lognorm, aic_lognorm2), 1)

# Pela comparação de aic, temos que o ajuste retirando as covariáveis que não
# apresentam significância sendo o melhor ajuste. Dessa forma seguiremos para a
# análise de resíduos
par(mfrow = c(1, 1))

plot(Ajuste_lognormal2)
# Apresenta boa aderência dos resíduo quantilicos

# a seed plantada no inicio do código foi necessária, pois nessa etapa o warm
# plot estava apresentando variação cada vez que era gerado, nesta seed o 
# gráfico apresentado no artigo foi constante
wp(Ajuste_lognormal2)

# Todas as observações estão dentro das bandas de confiança, indicando que o 
# ajuste permanece adequado

# Vamos obter o extrato1
exp(coef(Ajuste_lognormal2))[-1]

# Vamos obter o extrato2
exp(-coef(Ajuste_lognormal2))[-1]

# Curiosidade: o ajuste que desconsidera o efeito do tamanho do fígado acabou sendo mais
# útil, então na presença das demais variáveis o tamanho do fígado não aparenta
# conseguir contribuir na resposta.

# Realizei um ajuste considerando a variável do tamanho do fígado, porém não se
# mostrou melhor que o ajuste anterior.


