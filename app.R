# =========================================================================
# SOYRISK PREDICTOR: INTELIGÊNCIA ATUARIAL E FITOSSANITÁRIA
# =========================================================================

# -------------------------------------------------------------------------
# BLOCO 1: IMPORTAÇÃO DE PACOTES E DEPENDÊNCIAS
# -------------------------------------------------------------------------
library(shiny)          # Framework base para criar o aplicativo web interativo.
library(dplyr)          # Manipulação e transformação de dados (verbos como mutate, filter).
library(nasapower)      # Conexão com a API da NASA para extrair dados climáticos diários.
library(lubridate)      # Manipulação avançada de datas (cálculo de ciclo fenológico).
library(zoo)            # Aplicação de funções em janelas móveis (rolling windows) em séries temporais.
library(bslib)          # Estilização da interface de usuário utilizando Bootstrap 5.
library(leaflet)        # Renderização de mapas interativos profissionais.
library(leaflet.extras) # Extensão do leaflet para inclusão de Heatmaps (mapas de calor).
library(geobr)          # Acesso às malhas espaciais oficiais do IBGE mantidas pelo IPEA.
library(sf)             # (Simple Features) Processamento matemático de geometrias e polígonos espaciais.

# -------------------------------------------------------------------------
# BLOCO 2: INGESTÃO DE DADOS ESPACIAIS E MECANISMO SAFE-LOAD
# -------------------------------------------------------------------------
usando_fallback <- FALSE

# 2.1. Extração da malha estadual
estados_br_raw <- tryCatch({
  suppressMessages(geobr::read_state(year = 2020, showProgress = FALSE))
}, error = function(e) NULL)

if (!is.null(estados_br_raw) && nrow(estados_br_raw) > 0) {
  estados_br <- estados_br_raw %>% sf::st_transform(crs = 4326)
  estados_vizinhos <- estados_br %>% filter(abbrev_state != "GO")
  borda_goias <- estados_br %>% filter(abbrev_state == "GO")
} else {
  estados_vizinhos <- sf::st_sf(geometry = sf::st_sfc(), crs = 4326)
  borda_goias <- sf::st_sf(geometry = sf::st_sfc(), crs = 4326)
}

# 2.2. Extração da malha municipal (Estado de Goiás)
municipios_geo_raw <- tryCatch({
  suppressMessages(geobr::read_municipality(code_muni = "GO", year = 2020, showProgress = FALSE))
}, error = function(e) NULL)

if (!is.null(municipios_geo_raw) && nrow(municipios_geo_raw) > 0) {
  municipios_geo <- municipios_geo_raw %>% sf::st_transform(crs = 4326)
  centroides <- suppressWarnings(sf::st_coordinates(sf::st_centroid(municipios_geo)))
  municipios_geo$longitude <- centroides[,1]
  municipios_geo$latitude <- centroides[,2]
} else {
  usando_fallback <- TRUE
  demo_cities <- data.frame(
    name_muni = c("Rio Verde", "Jataí", "Cristalina", "Catalão", "Formosa", "Porangatu", "Ipameri", "Mineiros"),
    latitude = c(-17.79, -17.88, -16.76, -18.16, -15.53, -13.44, -17.72, -17.56),
    longitude = c(-50.92, -51.71, -47.61, -47.94, -47.33, -48.22, -48.15, -52.55)
  )
  municipios_geo <- sf::st_as_sf(demo_cities, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
}

# 2.3. Simulação Estocástica de Base Produtiva (Geração de cenário)
set.seed(42)
municipios_geo <- municipios_geo %>%
  mutate(
    severidade_historica = case_when(
      longitude < -50 & latitude < -17 ~ runif(n(), 85, 100), 
      longitude > -48 ~ runif(n(), 40, 75),                   
      TRUE ~ runif(n(), 60, 90)                               
    ),
    area_soja_ha = runif(n(), 5000, 90000),
    produtividade_atingivel = runif(n(), 3200, 4100)
  )

municipios_choices <- sort(unique(municipios_geo$name_muni))
cores_escala <- c("#ffeb3b", "#ff9800", "#e53935", "#b71c1c", "#4a0007")

# -------------------------------------------------------------------------
# BLOCO 3: MOTOR DE CÁLCULO ATUARIAL E EPIDEMIOLÓGICO
# -------------------------------------------------------------------------
calcular_risco_app <- function(municipio_nome, lat, lon, data_plantio, area_soja_ha, gmr, nivel_cobertura, preco_saca) {
  
  data_plantio <- as.Date(data_plantio)
  data_R1 <- data_plantio + days(50) 
  dias_exposicao <- 35 + ((gmr - 6.0) * 16)
  data_fim_avaliacao <- data_R1 + days(round(dias_exposicao))
  
  clima_realtime <- tryCatch({
    get_power(
      community = "ag", lonlat = c(lon, lat), pars = c("PRECTOTCORR"),
      dates = c(as.character(data_R1), as.character(data_fim_avaliacao)), temporal_api = "daily"
    )
  }, error = function(e) stop("Erro de ligação com a API da NASA POWER."))
  
  janela_critica <- clima_realtime %>%
    arrange(YYYYMMDD) %>%
    mutate(
      chuva_30d = rollapplyr(PRECTOTCORR, 30, sum, fill = NA),
      dias_chuva_30d = rollapplyr(PRECTOTCORR > 0, 30, sum, fill = NA)
    ) %>%
    filter(!is.na(chuva_30d)) 
  
  if (nrow(janela_critica) == 0) {
    rain_acc <- sum(clima_realtime$PRECTOTCORR); rain_nod <- sum(clima_realtime$PRECTOTCORR > 0)
  } else {
    janela_critica <- janela_critica %>% arrange(desc(dias_chuva_30d), desc(chuva_30d)) %>% slice(1)
    rain_acc <- janela_critica$chuva_30d; rain_nod <- janela_critica$dias_chuva_30d
  }
  
  rain_acc_limitada <- min(rain_acc, 600)
  sev_estimada <- -3.8983 + (0.3777 * rain_acc_limitada) - (0.0003 * (rain_acc_limitada^2))
  sev_estimada <- min(max(sev_estimada, 0), 100) 
  
  prod_esperada_kgha <- 3600; coeficiente_dano_kgha <- 21.41
  perda_kg_ha <- min(sev_estimada * coeficiente_dano_kgha, prod_esperada_kgha)
  prod_real_kgha <- prod_esperada_kgha - perda_kg_ha
  perda_total_toneladas <- (perda_kg_ha * area_soja_ha) / 1000
  
  prod_garantida_kgha <- prod_esperada_kgha * (nivel_cobertura / 100)
  if (prod_real_kgha < prod_garantida_kgha) {
    indenizacao_kg_ha <- prod_garantida_kgha - prod_real_kgha
    status_seguro <- "SINISTRO CONFIRMADO"; cor_status <- "#d9534f"
  } else {
    indenizacao_kg_ha <- 0; status_seguro <- "SEM INDENIZAÇÃO"; cor_status <- "#5cb85c"
  }
  
  custo_total_seguradora_brl <- (indenizacao_kg_ha / 60) * preco_saca * area_soja_ha 
  
  tibble(
    Municipio = municipio_nome, Lat = lat, Lon = lon, GMR = gmr, Cobertura_Pct = nivel_cobertura, Preco_Saca_Simulado = preco_saca,
    Data_Plantio = data_plantio, Data_R1 = data_R1, Data_Fim = data_fim_avaliacao, Dias_Exposicao = round(dias_exposicao),
    Chuva_Acumulada_30d_mm = round(rain_acc, 1), Dias_Chuvosos_30d = rain_nod, Severidade_Ferrugem_Predita = round(sev_estimada, 2),
    Perda_Estimada_kg_ha = round(perda_kg_ha, 2), Prejuizo_Fisico_Sc_Ha = round(perda_kg_ha / 60, 1),
    Prod_Real_Est_kg_ha = round(prod_real_kgha, 2), Perda_Total_Municipio_Tons = round(perda_total_toneladas, 2),
    Status_Seguro = status_seguro, Cor_Status = cor_status, Indenizacao_Total_Estimada_BRL = round(custo_total_seguradora_brl, 2),
    Nivel_de_Risco_Biologico = case_when(sev_estimada >= 75 ~ "CRÍTICO", sev_estimada >= 45 ~ "ALTO", sev_estimada >= 20 ~ "MÉDIO", TRUE ~ "BAIXO")
  )
}

# -------------------------------------------------------------------------
# BLOCO 4: INTERFACE DE USUÁRIO (UI)
# -------------------------------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "litera", primary = "#2E7D32", base_font = font_google("Montserrat")),
  
  tags$head(
    tags$style(HTML("
      body, p, div, span, table, th, td, h1, h2, h3, h4, h5, h6 { font-family: 'Montserrat', sans-serif !important; }
      .well { background-color: #f8f9fa; border: none; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); }
      .btn-primary { border-radius: 8px; text-transform: uppercase; letter-spacing: 1px; font-weight: bold; }
      .relatorio-card { padding: 25px; border-radius: 12px; background: white; box-shadow: 0 10px 20px rgba(0,0,0,0.08); border-top: 4px solid #2E7D32; margin-top: 20px; }
      .risco-box { padding: 15px; border-radius: 8px; color: white; font-weight: 700; text-align: center; font-size: 18px; margin: 20px 0; }
      .negocio-box { border: 2px solid #e9ecef; padding: 20px; border-radius: 8px; background-color: #ffffff; }
      .table-container { overflow-x: auto; background: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
      .table-striped tbody tr:nth-of-type(odd) { background-color: #f6f8fa; }
      .radio-btn-group .radio { padding: 10px; border: 1px solid #ddd; border-radius: 8px; background: white; margin-bottom: 5px; transition: 0.3s; }
      .radio-btn-group .radio:hover { background: #f1f8e9; border-color: #2E7D32; }
    "))
  ),
  
  titlePanel(
    div(icon("leaf", class = "text-primary"), "SoyRisk Predictor - Goiás",
        span("Mapeamento e Monitoramento Atuarial Avançado", style="font-size: 16px; color: #7f8c8d; margin-left: 15px;"))
  ),
  
  br(),
  
  sidebarLayout(
    sidebarPanel(
      h4(icon("power-off"), " Modo de Avaliação"),
      div(class = "radio-btn-group",
          radioButtons("modo_analise", label = NULL, 
                       choices = c("Completo (Risco Biológico + Financeiro)" = "completo", 
                                   "Apenas Risco Biológico (Fitopatologia)" = "bio"), 
                       selected = "completo")
      ),
      hr(),
      
      h4(icon("sliders-h"), " Parâmetros do Talhão"),
      selectInput("municipio", "Município Alvo", choices = municipios_choices, selected = "Rio Verde"),
      
      div(style = "display: none;",
          numericInput("lat", "Latitude", value = -17.79, step = 0.01),
          numericInput("lon", "Longitude", value = -50.92, step = 0.01)
      ),
      
      dateInput("data_plantio", "Data de Plantio", value = "2024-10-20"),
      numericInput("area_soja", "Área Plantada (Hectares)", value = 5000, step = 100),
      sliderInput("gmr", "Grau de Maturação Relativa (GMR)", min = 6.0, max = 8.5, value = 7.5, step = 0.1),
      
      conditionalPanel(
        condition = "input.modo_analise == 'completo'",
        hr(),
        h4(icon("shield-alt"), " Parâmetros da Seguradora"),
        sliderInput("cobertura", "Nível de Cobertura da Apólice (%)", min = 60, max = 85, value = 70, step = 5),
        sliderInput("preco_saca", "Preço Simulado da Saca", min = 90, max = 200, value = 135, step = 1, pre = "R$ ")
      ),
      
      br(),
      actionButton("calcular", "Executar Análise", class = "btn-primary w-100 p-3")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Dashboard de Risco",
                 br(),
                 h4(icon("map-marked-alt"), " Painel Cartográfico de Riscos Biológicos"),
                 leafletOutput("heatmap_estado", height = "520px"),
                 uiOutput("relatorio_ui"),
                 
                 br(),
                 h4(icon("list-alt"), " Matriz de Auditoria (Justificativa de Cálculo)"),
                 hr(),
                 div(class = "table-container", tableOutput("tabela_resultados")),
                 
                 br(), br(),
                 h4(icon("chart-bar"), " Portfólio de Risco: Top 10 Municípios de Goiás"),
                 hr(),
                 div(class = "table-container", tableOutput("tabela_ranking"))
        ),
        
        tabPanel("Metodologia & Código",
                 br(),
                 h4(icon("book"), " Fundamentação Científica"),
                 hr(),
                 tags$ul(
                   tags$li(tags$b("Sistema de Alerta Baseado em Chuva: "), "Beruski, G. C., Del Ponte, E. M., et al. (2020). Performance and profitability of rain-based thresholds for timing fungicide applications in soybean rust control. ", tags$em("Plant Disease, 104"), "(10), 2704-2712."),
                   tags$li(tags$b("Modelagem Empírica (Equação BR3): "), "Del Ponte, E. M., Godoy, C. V., Li, X., & Yang, X. B. (2006). Predicting severity of Asian soybean rust epidemics with empirical rainfall models. ", tags$em("Phytopathology, 96"), "(7), 797-803."),
                   tags$li(tags$b("Perda Financeira vs. Severidade: "), "Dalla Lana, F., Ziegelmann, P. K., de Maia, A. H. N., Godoy, C. V., & Del Ponte, E. M. (2015). Meta-analysis of the relationship between crop yield and soybean rust severity. ", tags$em("Phytopathology, 105"), "(3), 307-315."),
                   tags$li(tags$b("Ancoragem Fenológica (Escala Estágios Reprodutivos): "), "Fehr, W. R., & Caviness, C. E. (1977). Stages of soybean development. ", tags$em("Special Report 80"), ". Iowa State University.")
                 ),
                 br(),
                 h4(icon("book-open"), " Parâmetros e Normativas Regionais"),
                 hr(),
                 tags$ul(
                   tags$li(tags$b("Embrapa Soja e Parceiros: "), "EMBRAPA. Tecnologias de Produção de Soja - Região Central do Brasil."),
                   tags$li(tags$b("Reuniões Oficiais do Cerrado: "), "EMBRAPA SOJA. Atas da Reunião de Pesquisa de Soja da Região Central do Brasil."),
                   tags$li(tags$b("Universidade Federal de Goiás (UFG): "), "Pesquisa Agropecuária Tropical (ISSN: 1983-4063)."),
                   tags$li(tags$b("Zoneamento Agrícola: "), "MAPA / Embrapa. Zoneamento Agrícola de Risco Climático (ZARC).")
                 ),
                 br(),
                 h4(icon("users"), " Sobre o Time"),
                 hr(),
                 tags$ul(
                   tags$li(tags$b("Enzo Agustin Pedraza | "), tags$a(href = "https://www.linkedin.com/in/enzo-agustin-pedraza-6ab8b734a/", target = "_blank", "LinkedIn"), tags$br(), "Engenheiro Agrônomo formado pela Universidad Nacional de Tucumán (Argentina). Atualmente, é integrante do Laboratório de Epidemiologia da UFV."),
                   tags$li(tags$b("Gabriel José Degaspari | "), tags$a(href = "https://www.linkedin.com/in/gabriel-jos%C3%A9-degaspari-2720a5263/", target = "_blank", "LinkedIn"), tags$br(), "Engenheiro Agrônomo graduado pela própria Universidade Federal de Viçosa (Brasil). Desenvolve seu projeto de mestrado no Laboratório de Nematologia da UFV."),
                   tags$li(tags$b("Laura Gomez Agudelo | "), tags$a(href = "https://www.linkedin.com/in/laura-gomez-ag27/", target = "_blank", "LinkedIn"), tags$br(), "Engenheira Agropecuária graduada pelo Politécnico Colombiano 'Jaime Isaza Cadavid' (Colômbia). Integra o Laboratório de Controle Biológico de Nematoides da UFV."),
                   tags$li(tags$b("Milton Epitácio Carneiro Monte Galvino | "), tags$a(href = "https://www.linkedin.com/in/milton-monte-galvino-345674155/", target = "_blank", "LinkedIn"), tags$br(), "Engenheiro Agrônomo formado pela Universidade Federal do Ceará (Brasil). Atua como integrante do Laboratório de Biologia de Populações de Fitopatógenos da UFV.")
                 ),
                 br(),
                 h4(icon("github"), " Repositório do Projeto (GitHub)"),
                 hr(),
                 tags$a(href = "https://github.com/miltongalvino/SoyRisk-Predictor---Goi-s.git", target = "_blank", style = "font-size: 18px; font-weight: bold; color: #2E7D32;", "github.com/miltongalvino")
        )
      )
    )
  )
)

# -------------------------------------------------------------------------
# BLOCO 5: LÓGICA DE SERVIDOR (SERVER)
# -------------------------------------------------------------------------
server <- function(input, output, session) {
  
  observeEvent(input$municipio, {
    muni_data <- subset(municipios_geo, name_muni == input$municipio)
    if (nrow(muni_data) > 0) {
      updateNumericInput(session, "lat", value = muni_data$latitude[1])
      updateNumericInput(session, "lon", value = muni_data$longitude[1])
    }
  })
  
  output$heatmap_estado <- renderLeaflet({
    pal_legenda <- colorNumeric(palette = cores_escala, domain = c(40, 100))
    
    mapa <- leaflet(options = leafletOptions(minZoom = 6.0, maxZoom = 13)) %>%
      addProviderTiles(providers$CartoDB.Positron, group = "Mapa Claro (Foco)") %>% 
      addProviderTiles(providers$Esri.WorldImagery, group = "Visão de Satélite") %>%
      addProviderTiles(providers$OpenStreetMap, group = "Mapa Técnico (OSM)") %>%
      setView(lng = -49.6, lat = -16.0, zoom = 6.8) %>%
      setMaxBounds(lng1 = -54.0, lat1 = -20.0, lng2 = -45.0, lat2 = -12.0) %>%
      addHeatmap(
        data = municipios_geo, 
        lng = ~longitude, 
        lat = ~latitude, 
        intensity = ~severidade_historica, 
        blur = 12, 
        max = 60, 
        radius = 25, 
        gradient = cores_escala, 
        group = "Mancha de Severidade (Heatmap)"
      )
    
    if (!usando_fallback) {
      mapa <- mapa %>% 
        addPolygons(
          data = municipios_geo, 
          fillColor = ~pal_legenda(severidade_historica), 
          fillOpacity = 0.6, 
          color = "#ffffff", 
          weight = 0.7, 
          stroke = TRUE, 
          group = "Severidade por Município", 
          popup = ~paste0("<div style='font-family: Montserrat; padding: 5px;'><h6 style='font-weight:bold; color:#2E7D32;'>", name_muni, "</h6><b>Severidade Histórica:</b> ", round(severidade_historica, 1), "%</div>")
        ) %>%
        addPolygons(
          data = estados_vizinhos, 
          fillColor = "#e6e6e6", 
          fillOpacity = 0.85, 
          color = "#b0b0b0", 
          weight = 1.0, 
          stroke = TRUE, 
          group = "Isolar Estado de Goiás"
        )
    } else {
      mapa <- mapa %>% 
        addCircleMarkers(
          data = municipios_geo, 
          lng = ~longitude, 
          lat = ~latitude, 
          radius = 12, 
          fillColor = ~pal_legenda(severidade_historica), 
          color = "#fff", 
          weight = 1, 
          fillOpacity = 0.8, 
          group = "Severidade por Município", 
          popup = ~paste0("<div style='font-family: Montserrat; padding: 5px;'><h6 style='font-weight:bold; color:#2E7D32;'>", name_muni, "</h6><b>Severidade Histórica:</b> ", round(severidade_historica, 1), "%</div>")
        )
    }
    
    mapa <- mapa %>% 
      addPolygons(
        data = borda_goias, 
        fillColor = "transparent", 
        color = "#2E7D32", 
        weight = 2, 
        opacity = 0.9, 
        stroke = TRUE, 
        group = "Limite Estadual (GO)"
      ) %>%
      addLegend(
        pal = pal_legenda, 
        values = c(40, 100), 
        title = "Severidade Histórica (%)", 
        position = "bottomright", 
        opacity = 0.8
      ) %>%
      addMiniMap(
        tiles = providers$CartoDB.Positron, 
        toggleDisplay = TRUE, 
        position = "bottomleft"
      ) %>%
      addLayersControl( 
        baseGroups = c("Mapa Claro (Foco)", "Visão de Satélite", "Mapa Técnico (OSM)"),
        overlayGroups = c("Mancha de Severidade (Heatmap)", "Severidade por Município", "Isolar Estado de Goiás", "Limite Estadual (GO)"),
        options = layersControlOptions(collapsed = FALSE)
      ) %>%
      hideGroup("Severidade por Município")
    
    return(mapa)
  })
  
  dados_risco <- eventReactive(input$calcular, {
    withProgress(message = 'Consultando NASA POWER...', value = 0, {
      calcular_risco_app(municipio_nome = input$municipio, lat = input$lat, lon = input$lon, data_plantio = input$data_plantio, area_soja_ha = input$area_soja, gmr = as.numeric(input$gmr), nivel_cobertura = as.numeric(input$cobertura), preco_saca = as.numeric(input$preco_saca))
    })
  })
  
  observeEvent(input$calcular, {
    leafletProxy("heatmap_estado") %>% clearMarkers() %>%
      addCircleMarkers(lng = input$lon, lat = input$lat, radius = 9, color = "#1a237e", fillColor = "#00e5ff", fillOpacity = 1, weight = 2.5) %>%
      flyTo(lng = input$lon, lat = input$lat, zoom = 9)
  })
  
  output$tabela_resultados <- renderTable({
    req(dados_risco())
    df <- dados_risco()
    
    df_vertical <- tibble(
      `Parâmetro / Métrica` = c("Município Alvo", "Grau de Maturação Relativa (GMR)", "Nível de Cobertura Contratual", "Data de Plantio", "Data de R1 (Início de Risco)", "Dias Sob Risco (Exposição)", "Término da Janela Crítica", "Chuva Acumulada (30d)", "Dias Chuvosos (30d)", "Severidade Predita da Ferrugem", "Perda Física Estimada por Hectare", "Produtividade Real Projetada", "Status de Sinistro Atuarial", "Indemnização Estimada Total"),
      `Valor Estimado` = c(df$Municipio, as.character(df$GMR), paste0(df$Cobertura_Pct, "%"), format(df$Data_Plantio, "%d/%m/%Y"), format(df$Data_R1, "%d/%m/%Y"), paste(df$Dias_Exposicao, "dias"), format(df$Data_Fim, "%d/%m/%Y"), paste(df$Chuva_Acumulada_30d_mm, "mm"), paste(df$Dias_Chuvosos_30d, "dias"), paste0(df$Severidade_Ferrugem_Predita, "%"), paste(format(df$Perda_Estimada_kg_ha, big.mark = ".", decimal.mark = ","), "kg/ha"), paste(format(df$Prod_Real_Est_kg_ha, big.mark = ".", decimal.mark = ","), "kg/ha"), df$Status_Seguro, paste("R$", format(df$Indenizacao_Total_Estimada_BRL, big.mark = ".", decimal.mark = ","))),
      `Metodologia de Cálculo` = c("Geocodificação IBGE", "Características da Cultivar", "Garantia PROAGRO/PSR", "Data de Semeadura", "Data de Plantio + 50 dias", "35 dias + ((GMR - 6.0) x 16)", "Data de R1 + Dias de Exposição ao Risco", "Soma pluviométrica (NASA)", "Contagem de dias c/ chuva > 0mm", "-3.8983 + (0.3777 x Chuva_mm) - (0.0003 x Chuva_mm²)", "Severidade % x Coeficiente de Dano (21.41)", "Produtividade Base Esperada (3600 kg/ha) - Perda Física", "Se Prod. Real < (Prod. Base x Cobertura %)", "(Prod. Garantida - Prod. Real) x Área x Preço / 60"),
      `Fonte Bibliográfica` = c("Pacote geobr / IBGE", "Critério do Produtor", "Apólice / Seguradora", "Vazio Sanitário (Agrodefesa)", "Fehr & Caviness (1977)", "EMBRAPA Soja (Sist. de Produção)", "UFG / Pesq. Agrop. Trop.", "API NASA POWER", "API NASA POWER", "Modelo BR3 (Del Ponte et al., 2006)", "Dalla Lana et al. (2015)", "Cálculo Interno", "Regras SUSEP / MAPA", "Regras SUSEP / MAPA")
    )
    
    if (input$modo_analise == "bio") {
      df_vertical <- df_vertical %>%
        filter(!(`Parâmetro / Métrica` %in% c("Nível de Cobertura Contratual", "Produtividade Real Projetada", "Status de Sinistro Atuarial", "Indemnização Estimada Total")))
    }
    
    return(df_vertical)
  }, sanitize.text.function = function(x) x)
  
  output$tabela_ranking <- renderTable({
    req(input$preco_saca)
    
    df_ranking <- municipios_geo %>% sf::st_drop_geometry() %>%
      mutate(
        Perda_kg_ha = ifelse(severidade_historica * 21.41 > produtividade_atingivel, produtividade_atingivel, severidade_historica * 21.41),
        Perda_Total_Ton = (Perda_kg_ha * area_soja_ha) / 1000,
        Prejuizo_Financeiro_BRL = (Perda_Total_Ton * 1000 / 60) * input$preco_saca
      )
    
    if (input$modo_analise == "bio") {
      df_ranking %>% arrange(desc(Perda_Total_Ton)) %>% head(10) %>%
        mutate(
          severidade_historica = paste0(round(severidade_historica, 1), "%"), 
          area_soja_ha = format(round(area_soja_ha, 0), big.mark = ".", decimal.mark = ","), 
          produtividade_atingivel = format(round(produtividade_atingivel, 0), big.mark = ".", decimal.mark = ","), 
          Perda_Total_Ton = format(round(Perda_Total_Ton, 0), big.mark = ".", decimal.mark = ",")
        ) %>%
        select(name_muni, severidade_historica, area_soja_ha, produtividade_atingivel, Perda_Total_Ton) %>%
        rename(
          `Município` = name_muni, 
          `Severidade Biológica` = severidade_historica, 
          `Superfície Simulada (ha)` = area_soja_ha, 
          `Produtividade Teto (kg/ha)` = produtividade_atingivel, 
          `Perda Física Estimada (Toneladas)` = Perda_Total_Ton
        )
    } else {
      df_ranking %>% arrange(desc(Prejuizo_Financeiro_BRL)) %>% head(10) %>%
        mutate(
          severidade_historica = paste0(round(severidade_historica, 1), "%"), 
          area_soja_ha = format(round(area_soja_ha, 0), big.mark = ".", decimal.mark = ","), 
          produtividade_atingivel = format(round(produtividade_atingivel, 0), big.mark = ".", decimal.mark = ","), 
          Perda_Total_Ton = format(round(Perda_Total_Ton, 0), big.mark = ".", decimal.mark = ","), 
          Prejuizo_Financeiro_BRL = paste("R$", format(round(Prejuizo_Financeiro_BRL, 2), big.mark = ".", decimal.mark = ","))
        ) %>%
        select(name_muni, severidade_historica, area_soja_ha, produtividade_atingivel, Perda_Total_Ton, Prejuizo_Financeiro_BRL) %>%
        rename(
          `Município` = name_muni, 
          `Severidade Biológica` = severidade_historica, 
          `Superfície Simulada (ha)` = area_soja_ha, 
          `Produtividade Teto (kg/ha)` = produtividade_atingivel, 
          `Perda Física (Toneladas)` = Perda_Total_Ton, 
          `Risco Financeiro Estimado` = Prejuizo_Financeiro_BRL
        )
    }
  })
  
  output$relatorio_ui <- renderUI({
    req(dados_risco())
    res <- dados_risco()
    
    cor_risco <- switch(res$Nivel_de_Risco_Biologico, "CRÍTICO" = "#d9534f", "ALTO" = "#f0ad4e", "MÉDIO" = "#f39c12", "BAIXO" = "#2ecc71")   
    
    html_biologico <- sprintf(
      "<div class='relatorio-card'>
        <h5><i class='fa fa-microscope text-primary'></i> Relatório de Avaliação Fenológica e Climática</h5>
        <hr>
        <p><b>Município Enquadrado:</b> %s, GO | <b>Variedade:</b> GMR %s (%s dias sob risco)</p>
        <p><b>Janela Alvo de Monitoramento:</b> %s até %s</p>
        <div class='risco-box' style='background-color: %s;'>
          Risco de Perda: %s (Severidade Real Medida: %s%%)
        </div>",
      res$Municipio, res$GMR, res$Dias_Exposicao, format(res$Data_R1, "%d/%m/%Y"), format(res$Data_Fim, "%d/%m/%Y"),
      cor_risco, res$Nivel_de_Risco_Biologico, res$Severidade_Ferrugem_Predita
    )
    
    if (input$modo_analise == "completo") {
      html_financeiro <- sprintf(
        "<div class='negocio-box' style='border-color: %s;'>
          <h5 style='color: %s; font-weight:700;'><i class='fa fa-dollar-sign'></i> Impacto Comercial Atuarial</h5>
          <p>Considerando o valor de mercado simulado a <b>R$ %s</b> por saca de 60kg e Cobertura de <b>%s%%</b>.</p>
          <p><b>Produtividade Estimada pós-infecção:</b> %s kg/ha (Perda física de %s sc/ha)</p>
          <div style='background: #f8f9fa; padding: 12px; border-radius: 6px; text-align: center; margin-top: 10px;'>
            <span style='font-size: 14px; color: #555;'>Previsão de Indemnização Financeira:</span>
            <h4 style='color: %s; margin: 5px 0 0 0; font-weight: 700;'>R$ %s</h4>
          </div>
        </div>
        </div>",
        res$Cor_Status, res$Cor_Status, format(res$Preco_Saca_Simulado, big.mark = ".", decimal.mark = ","), res$Cobertura_Pct,
        format(res$Prod_Real_Est_kg_ha, big.mark = ".", decimal.mark = ","), res$Prejuizo_Fisico_Sc_Ha,
        res$Cor_Status, format(res$Indenizacao_Total_Estimada_BRL, big.mark = ".", decimal.mark = ",")
      )
      HTML(paste0(html_biologico, html_financeiro))
    } else {
      HTML(paste0(html_biologico, "</div>"))
    }
  })
}

shinyApp(ui = ui, server = server)