
# #  SoyRisk Predictor - Goiás

### *Inteligência Atuarial e Fitossanitária Avançada para a Cultura da Soja*

O **SoyRisk Predictor** é uma aplicação web interativa desenvolvida em **R (Shiny)** que integra modelagem fitopatológica preditiva e análise de risco atuarial para a cultura da soja no estado de Goiás. O sistema monitora a janela crítica do ciclo fenológico da cultura e consome dados climáticos em tempo real para estimar a severidade da **Ferrugem Asiática da Soja** (*Phakopsora pachyrhizi*), traduzindo o risco biológico em impacto financeiro e probabilidade de sinistro para apólices de seguro agrícola (PROAGRO/PSR).

---

##  Funcionalidades Principais

* **Ingestão de Dados em Tempo Real:** Conexão automatizada com a API **NASA POWER** para coleta de dados pluviométricos diários com base nas coordenadas geográficas do município selecionado.
* **Análise Fenológica Personalizada:** Ajuste dinâmico do período de exposição ao risco a partir do Grau de Maturação Relativa (GMR) da cultivar utilizada.
* **Motor Epidemiológico Integrado:** Simulação da severidade da Ferrugem Asiática usando modelos empíricos consagrados na literatura científica.
* **Painel Cartográfico Interativo:** Visualização espacial com mapas de calor (Heatmaps) e polígonos municipais via `leaflet` e malhas oficiais do IBGE (`geobr`).
* **Simulação Atuarial Dinâmica:** Estimativa de perda física (sc/ha), projeção de produtividade real pós-infecção, validação de quebra de franquia/cobertura e cálculo de indenização total em moeda corrente (R$).
* **Dual Mode de Visualização:** Permite alternar entre uma análise completa (Biológica + Financeira) ou focar estritamente nos parâmetros fitopatológicos.

---

##  Tecnologias e Dependências

O projeto foi construído utilizando o ecossistema R, dependendo dos seguintes pacotes:

| Pacote | Função no Projeto |
| --- | --- |
| `shiny` | Framework base para a interface web interativa. |
| `bslib` | Estilização moderna da UI com Bootstrap 5 (Tema *Litera* / Fonte *Montserrat*). |
| `dplyr` & `lubridate` | Manipulação de dados e tratamento avançado de janelas temporais e datas. |
| `nasapower` | Cliente de comunicação com a API de Dados Climatológicos da NASA. |
| `zoo` | Processamento de funções em janelas móveis (rolling windows) de 30 dias. |
| `leaflet` & `leaflet.extras` | Construção de mapas interativos e renderização de camadas de calor. |
| `geobr` & `sf` | Download de malhas espaciais do IBGE e processamento de geometrias vetorizadas. |

---

##  Fundamentação Teórica e Algoritmos

### 1. Janela de Exposição Crítica

O período em que a planta fica mais suscetível a danos severos inicia-se no estágio reprodutivo $R_1$ (estimado em $\text{Data de Plantio} + 50\text{ dias}$) e sua duração é modelada linearmente com base no ciclo da cultivar ($GMR$):

$$\text{Dias de Exposição} = 35 + ((GMR - 6.0) \times 16)$$

### 2. Modelo Epidemiológico (Modelo BR3)

A severidade da doença é estimada encontrando a janela móvel de 30 dias mais favorável (maior acúmulo e frequência de chuva) dentro do período de exposição. Aplica-se a equação empírica de *Del Ponte et al. (2006)*:

$$\text{Severidade (\%)} = -3.8983 + (0.3777 \times x) - (0.0003 \times x^2)$$

>  *Onde x representa a chuva acumulada da janela crítica limitada a 600 mm. A perda de produtividade física é calculada com base no coeficiente de dano médio de 21.41 kg/ha para cada 1% de severidade (Dalla Lana et al., 2015).*

---

##  Como Executar o Projeto

### Pré-requisitos

Certifique-se de ter o **R** e o **RStudio** instalados em sua máquina.

### Passos para Execução

1. Clone este repositório em sua máquina local:
```bash
git clone https://github.com/miltongalvino/SoyRisk-Predictor---Goias.git

```


2. Abra o arquivo do script no RStudio.
3. Instale as dependências necessárias executando o código abaixo no console do R:
```R
install.packages(c("shiny", "dplyr", "nasapower", "lubridate", "zoo", 
                   "bslib", "leaflet", "leaflet.extras", "geobr", "sf"))

```


4. Execute o aplicativo clicando no botão **"Run App"** no RStudio ou digitando:
```R
shiny::runApp()

```



---

##  Referências Científicas Utilizadas

* **Beruski, G. C., Del Ponte, E. M., et al. (2020).** Performance and profitability of rain-based thresholds for timing fungicide applications in soybean rust control. *Plant Disease*, 104(10), 2704-2712.
* **Del Ponte, E. M., Godoy, C. V., Li, X., & Yang, X. B. (2006).** Predicting severity of Asian soybean rust epidemics with empirical rainfall models. *Phytopathology*, 96(7), 797-803.
* **Dalla Lana, F., Ziegelmann, P. K., de Maia, A. H. N., Godoy, C. V., & Del Ponte, E. M. (2015).** Meta-analysis of the relationship between crop yield and soybean rust severity. *Phytopathology*, 105(3), 307-315.
* **Fehr, W. R., & Caviness, C. E. (1977).** Stages of soybean development. *Special Report 80*. Iowa State University.
* **MAPA / Embrapa.** Zoneamento Agrícola de Risco Climático (ZARC).

---

##  Equipe de Desenvolvimento (UFV)

Este projeto foi desenvolvido por engenheiros agrônomos e pesquisadores vinculados aos programas de pós-graduação da **Universidade Federal de Viçosa (UFV)**:

* **Enzo Agustin Pedraza** – Universidade Federal de Viçosa (BR). Integrante do Laboratório de Epidemiologia da UFV. [[LinkedIn](https://www.linkedin.com/in/enzo-agustin-pedraza-6ab8b734a/)]
* **Gabriel José Degaspari** – Universidade Federal de Viçosa (BR). Integrante do Laboratório de Nematologia da UFV. [[LinkedIn](https://www.linkedin.com/in/gabriel-jos%C3%A9-degaspari-2720a5263/)]
* **Laura Gomez Agudelo** – Universidade Federal de Viçosa (BR). Integra o Laboratório de Controle Biológico de Nematoides da UFV. [[LinkedIn](https://www.linkedin.com/in/laura-gomez-ag27/)]
* **Milton Epitácio Carneiro Monte Galvino** – Universidade Federal de Viçosa (BR). Atua no Laboratório de Biologia de Populações de Fitopatógenos da UFV. [[LinkedIn](https://www.linkedin.com/in/milton-monte-galvino-345674155/)]

---

🔗 **Repositório Oficial:** [github.com/miltongalvino/SoyRisk-Predictor---Goias](https://github.com/miltongalvino/SoyRisk-Predictor---Goias)
🔗 **Site do app**: https://miltongalvino.shinyapps.io/soyriskpredictor/
