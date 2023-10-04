library(shiny)
library(signal) # for Butterworth filter
library(ggplot2)
library(tidyverse) # for data manipulation
library(oce, quietly = T)
detect_emg_events <-
  function(emg_signal,
           sampling_rate,
           threshold = 0.1,
           window_size = 50,
           min_length = 0.5) {
    onset_times <- c()  # Initialize a vector to store onset times
    offset_times <- c()  # Initialize a vector to store offset times

    window_on <- 0   # Initialize a counter for the window size
    window_off <- 0
    in_onset <-
      FALSE   # Flag to keep track if we are already counting an onset
    marker <- -1 * (sampling_rate)
    for (i in 1:length(emg_signal)) {
      # If the EMG signal crosses the threshold (potential onset)
      if (emg_signal[i] > threshold & !in_onset) {
        window_on <- window_on + 1  # Increment the window counter

        # If the signal has stayed above the threshold for 'window_size' samples
        if (window_on == window_size &
            !in_onset &
            i >= (marker + sampling_rate * min_length)) {
          onset_time <-
            (i - window_size) / sampling_rate  # Calculate the onset time in seconds
          onset_times <-
            c(onset_times, onset_time)  # Append to the list of onset times
          in_onset <-
            TRUE  # Set the flag to indicate we are in an onset
          window_off <- 0
          marker <- i
        }

      } else if (emg_signal[i] < threshold &
                 in_onset &
                 i >= (marker + sampling_rate * min_length)) {
        # If the EMG signal goes below the threshold (potential offset)
        window_off <- window_off + 1  # Decrement the window counter

        # If the signal has stayed below the threshold for 'window_size' samples
        if (window_off == window_size & in_onset) {
          offset_time <-
            (i - window_size) / sampling_rate  # Calculate the offset time in seconds
          offset_times <-
            c(offset_times, offset_time)  # Append to the list of offset times
          in_onset <- FALSE  # Reset the onset flag
          window_on <- 0
          marker <- i
        }


      }
    }

    return(list(onset_times = onset_times, offset_times = offset_times))  # Return the onset and offset times
  }
# UI layout
ui <- fluidPage(titlePanel("EMG Signal Processing"),

                sidebarLayout(
                  sidebarPanel(
                    fileInput(
                      "file",
                      "Upload a CSV File",
                      accept = c("text/csv", "text/comma-separated-values, text/plain")
                    ),
                    numericInput("channel", "EMG channel:", 1, min = 1),
                    numericInput("samplingFreq", "Sampling Frequency (Hz):", 2000),
                    numericInput("filterOrder", "Bandpass Filter Order:", 4),
                    numericInput("lowpass", "Low pass filter frequency:", 20),
                    numericInput("highpass", "High pass filter frequency:", 450),
                    numericInput("envelopeFreq", "Envelope Filter Frequency (Hz):", 10),
                    checkboxInput("normalize", "Normalize EMG Amplitude?", FALSE),
                    conditionalPanel(
                      condition = "input.normalize == true",
                      textInput("normFactor", "Enter EMG Amplitude for Normalization:", "")
                    ),
                    sliderInput(
                      "timeRange",
                      "Time Frame for Analysis (s):",
                      min = 0,
                      max = 100,
                      value = c(0, 100),
                      step = 0.01,
                    ),
                    actionButton("update", "Update Plots")

                  ),

                  mainPanel(
                    tabsetPanel(
                      tabPanel(
                        "Time Domain",
                        actionButton("calc_iEMG", "Calculate EMG amplitudes"),
                        downloadButton("downloadData", "Download EMG Envelope"),
                        
                        plotOutput("plotOriginal"),
                        plotOutput("plotFiltered"),
                        plotOutput("plotEnvelope"),
                        verbatimTextOutput("iEMG")
                      ),
                      tabPanel(
                        "Frequency Domain",
                        downloadButton("downloadpower", "Download EMG power spectrum"),
                        actionButton("calc_medianFreq", "Calculate Median Frequencies"),
                        plotOutput("spectrogram"),
                        plotOutput("plotPowerSpectrum"),
                        verbatimTextOutput("overall_medianFreq"),
                        verbatimTextOutput("medianFreq")
                      ),
                      tabPanel(
                        "Automatic Segmentation",
                        numericInput("windowsize", "Window size:", 1000),
                        numericInput("threshold", "% peak EMG:", 5),
                        numericInput("min_gap", "Minimal duration of cycle (s):", 0.5),
                        plotOutput("on_off"),
                        actionButton("average1", "Block averaging: on - on"),
                        actionButton("average2", "Block averaging: on - off"),
                        plotOutput("average_plot"),
                        downloadButton("downloadonon", "Download on - on average"),
                        downloadButton("downloadonoff", "Download on - off average")
                      ),
                      tabPanel(
                        "Manual Segmentation",
                        textInput("event", "Time points of event markers:", '0,1,2'),
                        actionButton("segment", "Block averaging"),
                        plotOutput("segment_plot"),
                        plotOutput("average_cycle"),
                        downloadButton("downloadcycle", "Download cycle average")
                      ),
                      tabPanel(
                        "Latency",
                        numericInput("time_stim", "Time of Simulation:", 0),
                        actionButton("cal_delay", "Calculate Latency"),
                        plotOutput("delay_plot"),
                        verbatimTextOutput("latency")
                      )
                    )
                  )
                ))

# Server logic
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  inFile <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })
  envelopeSignal <- reactiveVal(NULL)
  medianFreqList <- reactiveVal(NULL)

  observe({
    if (is.null(inFile())) {
      return(NULL)
    }
    n <- input$channel + 1
    df <- data.frame(time = inFile()[, 1], signal = inFile()[, n])
    max_time <- max(inFile()[, 1], na.rm = TRUE)
    updateSliderInput(
      session,
      "timeRange",
      max = max_time,
      value = c(0, max_time),
      step = 0.01
    )

    updateNumericInput(session,
                       "channel",
                       max = ncol(inFile()) - 1,
                       min = 1)

    observe({
      selected_range <- input$timeRange
      start_time <- selected_range[1]
      end_time <- selected_range[2]
      output$plotOriginal <- renderPlot({
        ggplot(df, aes(x = time, y = signal)) +
          geom_line() +
          geom_vline(xintercept = start_time, col = 'red') +
          geom_vline(xintercept = end_time, col = 'blue') +
          ggtitle("Original Signal")
      })
      spec = specgram(
        x = df$signal,
        n = 512,
        Fs = input$samplingFreq,
        window = 256,
        overlap = 128
      )
      # discard phase information
      P = abs(spec$S)
      # normalize
      P = P / max(P)
      # convert to dB
      P = 10 * log10(P)
      # config time axis
      t = spec$t
      # plot spectrogram
      output$spectrogram <- renderPlot({
        imagep(
          x = t,
          y = spec$f,
          z = t(P),
          col = oce.colorsViridis,
          ylab = 'Frequency [Hz]',
          xlab = 'Time [s]',
          drawPalette = T,
          decimate = F,
          ylim = c(0, 500),
          main = 'Spectrogram'
        )
      })
    })

    observeEvent(input$update, {
      selected_range <- input$timeRange
      start_time <- selected_range[1]
      end_time <- selected_range[2]
      fs <- input$samplingFreq
      n <- input$filterOrder
      l <- input$lowpass
      h <- input$highpass
      Wn_bandpass <-
        c(l, h) / (fs / 2) # Normalize by Nyquist frequency
      selected_data <-
        df %>% dplyr::filter(time >= start_time & time <= end_time)
      dc_removed_signal <-
        selected_data$signal - mean(selected_data$signal)

      # Apply a band-pass Butterworth filter
      b_bandpass <- butter(n, Wn_bandpass, type = "pass")
      filtered_Signal <- filtfilt(b_bandpass, dc_removed_signal)

      output$plotFiltered <- renderPlot({
        ggplot(
          data.frame(
            time = 1:length(filtered_Signal) / fs,
            signal = filtered_Signal
          ),
          aes(x = time, y = signal)
        ) +
          geom_line() +
          ggtitle("Filtered Signal")
      })

      # Additional low-pass filter for EMG Envelope
      Wn_envelope <- input$envelopeFreq / (fs / 2)
      b_envelope <- butter(n, Wn_envelope, type = "low")
      envelopeSignal <-
        filtfilt(b_envelope, abs(filtered_Signal)) # rectify and filter
      envelopeSignal_copy <- envelopeSignal
      #EMG amplitude normalization. When the checkbox is checked
      if (input$normalize) {
        norm_factor <- as.numeric(input$normFactor)
        if (!is.na(norm_factor) && norm_factor != 0) {
          envelopeSignal <- envelopeSignal / norm_factor
        }
      }
      envelopeSignal(envelopeSignal)
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("EMG_envelope_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          data.frame(time = 1:length(envelopeSignal) / fs,
                     envelope = envelopeSignal) %>%
            write.csv(file = file, row.names = FALSE)
        }
      )
      output$plotEnvelope <- renderPlot({
        ggplot(
          data.frame(
            time = 1:length(filtered_Signal) / fs,
            signal = envelopeSignal
          ),
          aes(x = time, y = signal)
        ) +
          geom_line() +
          ggtitle("EMG Envelope")
      })
      ####event detection ###########

      onoff_event <-
        detect_emg_events(
          envelopeSignal,
          fs,
          window_size = input$windowsize,
          max(envelopeSignal) * input$threshold / 100,
          input$min_gap
        )
      if (!is.null(onoff_event[[1]]) & !is.null(onoff_event[[2]])) {
        output$on_off <- renderPlot({
          ggplot(
            data.frame(
              time = 1:length(filtered_Signal) / fs,
              signal = envelopeSignal
            ),
            aes(x = time, y = signal)
          ) +
            geom_line() + geom_vline(xintercept = onoff_event[[1]], col = 'red') +
            geom_vline(xintercept = onoff_event[[2]], col = 'blue') + ggtitle("EMG Envelope")
        })

      if (length(onoff_event[[1]]) > 1) {
        on <- onoff_event[[1]]
        n_cycle <- length(on) - 1
        cycle_mat <- matrix(NA, nrow = 1000, ncol = n_cycle)
        for (i in 1:n_cycle) {
          cycle <- envelopeSignal[(on[i] * fs):(on[i + 1] * fs)]
          cycle_mat[, i] <- approx(cycle, n  = 1000)$y
        }
        on_on_average <- rowMeans(as.data.frame(cycle_mat))
      }
      ######download on-on cycle######
      output$downloadonon <- downloadHandler(
        filename = function() {
          paste("Block average_cycle", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          data.frame(time = 1:1000 / 10,
                     EMG = on_on_average) %>%
            write.csv(file = file, row.names = FALSE)
        }
      )
      ####download on-off cycle#######

      on <- onoff_event[[1]]
      off <- onoff_event[[2]]
      n_off <- length(off)
      cycle_mat2 <- matrix(NA, nrow = 1000, ncol = n_off)
      for (i in 1:n_off) {
        cycle <- envelopeSignal[(on[i] * fs):(off[i] * fs)]
        cycle_mat2[, i] <- approx(cycle, n  = 1000)$y
      }
      on_off_average <- rowMeans(as.data.frame(cycle_mat2))

      output$downloadonoff <- downloadHandler(
        filename = function() {
          paste("Block average_onoff", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          data.frame(time = 1:1000 / 10,
                     EMG = on_off_average) %>%
            write.csv(file = file, row.names = FALSE)
        }
      )
      #####plot on-on block average##########
      observeEvent(input$average1, {
        if (length(onoff_event[[1]]) > 1) {
          output$average_plot <- renderPlot({
            ggplot(data.frame(#time = 1:length(cycle) / fs,
              time = (1:1000) / 10,
              signal = on_on_average),
              aes(x = time, y = signal)) + geom_line() +
              ggtitle("Average EMG Envelope") + xlab('% cycle')
          })
        }
      })
      #####plot on-off block average##########
      observeEvent(input$average2, {
        output$average_plot <- renderPlot({
          ggplot(data.frame(#time = 1:length(cycle) / fs,
            time = (1:1000) / 10,
            signal = on_off_average),
            aes(x = time, y = signal)) + geom_line() +

            ggtitle("Average EMG Envelope") + xlab('% cycle')
        })

      })
}
      ##########manual segmentation#####

      observeEvent(input$segment, {
        markers_location <-
          str_split(input$event, ',')[[1]] %>% as.numeric()
        output$segment_plot <- renderPlot({
          ggplot(
            data.frame(
              time = 1:length(filtered_Signal) / fs,
              signal = envelopeSignal
            ),
            aes(x = time, y = signal)
          ) +
            geom_line() + geom_vline(xintercept = markers_location, col = 'red') +
            ggtitle("EMG Envelope")
        })

        ncycle <- length(markers_location) - 1
        mcycle_mat <- matrix(NA, nrow = 1000, ncol = ncycle)
        for (i in 1:ncycle) {
          mcycle <-
            envelopeSignal[(markers_location[i] * fs):(markers_location[i + 1] * fs)]
          mcycle_mat[, i] <- approx(mcycle, n  = 1000)$y
        }
        mcycle_average <- rowMeans(as.data.frame(mcycle_mat))
        output$average_cycle <- renderPlot({
          ggplot(
            data.frame(time = 1:1000 / 10,
                       signal = mcycle_average),
            aes(x = time, y = signal)
          ) +
            geom_line() +
            xlab ("% cycle")
        })
        output$downloadcycle <- downloadHandler(
          filename = function() {
            paste("Block average_cycle", Sys.Date(), ".csv", sep = "")
          },
          content = function(file) {
            data.frame(time = 1:1000 / 10,
                       EMG =  mcycle_average) %>%
              write.csv(file = file, row.names = FALSE)
          }
        )
      })

      #####latency##########
      observeEvent(input$cal_delay, {
        on <- onoff_event[[1]][1]
        output$delay_plot <- renderPlot({
          ggplot(
            data.frame(
              time = 1:length(filtered_Signal) / fs,
              signal = envelopeSignal
            ),
            aes(x = time, y = signal)
          ) +
            geom_line() +
            ggtitle("EMG Envelope") + geom_point(aes(x = input$time_stim, y = 0),
                                                 col = 'red',
                                                 size = 5) +
            geom_point(aes(x = on, y = 0),
                       col = 'blue',
                       size = 5) +
            geom_point(aes(x = time[which.max(signal)], y = 0),
                       col = 'green',
                       size = 5)
        })
        time_to_peak <-
          which.max(envelopeSignal) / fs - input$time_stim
        output$latency <- renderText({
          paste(
            "Onset latency:",
            on[1] - input$time_stim,
            " s",
            '\n',
            'Stimulation to peak signal latency:',
            time_to_peak ,
            " s"
          )
        })
      })

      ##########calculate IEMG#######
      observeEvent(input$calc_iEMG, {
        iEMG_value <- sum(abs(envelopeSignal_copy)) * 1 / fs
        mean_EMG <- mean(abs(envelopeSignal_copy))
        peak_EMG <- max(abs(envelopeSignal_copy))


        output$iEMG <- renderText({
          paste(
            "Integrated EMG (iEMG) is:",
            round(iEMG_value, 7),
            '\n',
            "Mean EMG is:",
            mean_EMG,
            '\n',
            "Peak EMG is:",
            peak_EMG
          )
        })
      })
      observeEvent(input$calc_medianFreq, {
        fs <- input$samplingFreq
        selected_range <- input$timeRange
        start_time <- selected_range[1]
        end_time <- selected_range[2]
        one_sec_samples <- fs # number of samples in one second

        selected_data <-
          df %>% dplyr::filter(time >= start_time &
                                 time <= end_time)
        if (nrow(selected_data) == 0)
          return(NULL)

        dc_removed_signal <-
          filtfilt(b_bandpass,
                   selected_data$signal - mean(selected_data$signal))
        #dc_removed_signal <- inFile()$signal - mean(inFile()$signal)
        filtered_signal <- filtfilt(b_bandpass, dc_removed_signal)
        len <- length(filtered_signal)
        num_windows <- floor(len / one_sec_samples)

        median_frequencies <- numeric(num_windows)

        for (i in 1:num_windows) {
          start_idx <- (i - 1) * one_sec_samples + 1
          end_idx <- i * one_sec_samples
          window_data <- filtered_signal[start_idx:end_idx]

          pspectrum <- abs(fft(window_data)) ^ 2
          pspectrum <- pspectrum[1:(length(pspectrum) / 2)]

          total_power <- sum(pspectrum)
          cumulative_power <- cumsum(pspectrum)
          f <- seq(0, fs / 2, length.out = length(window_data) / 2)

          median_frequency <-
            f[which.min(abs(cumulative_power - total_power / 2))]
          median_frequencies[i] <- median_frequency
        }

        medianFreqList(median_frequencies)
      })

      # Display median frequencies
      output$medianFreq <- renderText({
        if (!is.null(medianFreqList())) {
          paste("Calculated Median Frequencies (Hz) on 1 second interval:",
                toString(round(medianFreqList(), 2)))
        } else {
          "Press the button to calculate median frequencies."
        }
      })

      ####### Power Spectrum Analysis#######
      f <- seq(0, fs / 2, length.out = length(filtered_Signal) / 2)
      pspectrum <- abs(fft(filtered_Signal)) ^ 2
      pspectrum <- pspectrum[1:(length(pspectrum) / 2)]
      length(f) <- length(pspectrum)

      output$plotPowerSpectrum <- renderPlot({
        ggplot(data.frame(frequency = f, spectrum = pspectrum),
               aes(x = frequency, y = spectrum)) +
          geom_line() +
          ggtitle("Power Spectrum")
      })

      output$downloadpower <- downloadHandler(
        filename = function() {
          paste("EMG_power_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          data.frame(frequency = f, spectrum = pspectrum) %>%
            write.csv(file = file, row.names = FALSE)
        }
      )

      # Median Frequency Calculation
      total_power <- sum(pspectrum)
      cumulative_power <- cumsum(pspectrum)
      median_frequency <-
        f[which.min(abs(cumulative_power - total_power / 2))]
      output$overall_medianFreq <- renderText({
        paste("Median Frequency is:",
              round(median_frequency, 3),
              "Hz")
      })

    })

  })

}

shinyApp(ui = ui, server = server)
