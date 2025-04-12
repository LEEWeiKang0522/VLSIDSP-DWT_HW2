# VLSIDSP-DWT_HW2
# Fixed-Point 2-D DWT Bit‑True Simulation

This project implements a fixed‑point simulation of the forward 2‑D Discrete Wavelet Transform (DWT) with a floating‑point Inverse DWT (IDWT) for image reconstruction. The goal is to determine the minimum word lengths needed so that the reconstructed image has a PSNR (Peak Signal-to-Noise Ratio) of at least 40 dB while using the smallest word lengths possible.

## Overview

The 2‑D DWT algorithm is performed in two steps:
1. **Horizontal 1‑D DWT:** Each row of an 8‑bit input image is processed by convolving with low‑pass and high‑pass filters and then downsampled.
2. **Vertical 1‑D DWT:** The outputs from the horizontal DWT are processed column‑wise to produce the final subbands (LL, LH, HL, HH).

For multi‑level processing (e.g., Levels 2 and 3), the LL subband from the previous level is used as input to the next level using the same fixed‑point format.

The fixed‑point forward DWT outputs are then converted to floating‑point using their defined scaling (LSB weight) before applying a floating‑point inverse DWT to reconstruct the image. PSNR is computed from the difference between the original and reconstructed images.

## Fixed‑Point Design Parameters

### A. Filter Coefficients
- **Word Length:** 16 bits  
- **Format:** 1‑bit sign + 1‑bit integer + 14‑bit fractional  
- **LSB Weight:** \(2^{-14}\)  

_All filter coefficients (low‑pass \(h(n)\) and high‑pass \(g(n)\)) are quantized to this uniform format._

### B. Intermediate Computations
- **Multiplication & Accumulation:**  
  Use 32‑bit accumulators to safely sum products from an 8‑bit input multiplied by a 16‑bit coefficient without overflow.

### C. DWT Output (Each Level)
- **Horizontal 1‑D DWT Output:**  
  - **Word Length:** 12 bits  
  - **Format:** 1‑bit sign + 3‑bit integer + 8‑bit fractional  
  - **LSB Weight:** \(2^{-8}\)  

- **Vertical 1‑D DWT (Final Subbands: LL, LH, HL, HH):**  
  - **Word Length:** 12 bits  
  - **Format:** Same as horizontal (1+3+8)  
  - **LSB Weight:** \(2^{-8}\)  

- **Multi‑Level DWT (Levels 2 & 3):**  
  Process the LL subband further using the same 12‑bit fixed‑point format.

## Bit‑True Algorithm Pseudocode

```c
// -------------------------------------------------
// Fixed-Point Parameter Definitions:
#define FILTER_WORD_LEN 16   // 16-bit coefficients: 1+1+14; LSB = 2^-14
#define ACCUM_WORD_LEN  32   // 32-bit accumulators for MAC operations
#define OUTPUT_WORD_LEN 12   // 12-bit DWT output: 1+3+8; LSB = 2^-8

// Pre-quantized filter coefficients (in 16-bit fixed point)
fixed16_t h[Nh] = { /* Example: quantize(0.852698679009, 16), ... */ };
fixed16_t g[Ng] = { /* Example: quantize(-0.788485616406, 16), ... */ };

// -------------------------------------------------
// Forward 2-D DWT: Horizontal 1-D DWT (Fixed-Point)
// -------------------------------------------------
for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j += step) {  // 'step' corresponds to downsampling factor
        int32_t acc_low  = 0;
        int32_t acc_high = 0;
        // Convolve over filter taps
        for (int k = -K; k <= K; k++) {
            fixed16_t sample = convert_to_fixed(input_image[i][j+k]);
            acc_low  += sample * h[k+offset];
            acc_high += sample * g[k+offset];
        }
        // Downsample and convert to 12-bit fixed point (with rounding/truncation)
        fixed12_t L[i][j/2] = truncate_round(acc_low, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
        fixed12_t H[i][j/2] = truncate_round(acc_high, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
    }
}

// -------------------------------------------------
// Forward 2-D DWT: Vertical 1-D DWT (Fixed-Point)
// -------------------------------------------------
for (int j = 0; j < (N/2); j++) {
    for (int i = 0; i < M; i += step) {
        int32_t acc_LL = 0;
        int32_t acc_LH = 0;
        // Convolve along the column
        for (int k = -K; k <= K; k++) {
            acc_LL += L[i+k][j] * h[k+offset];
            acc_LH += H[i+k][j] * g[k+offset];
        }
        // Downsample vertically and convert to 12-bit fixed point
        fixed12_t LL[i/2][j] = truncate_round(acc_LL, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
        fixed12_t LH[i/2][j] = truncate_round(acc_LH, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
        // Compute HL and HH similarly as needed
    }
}

// -------------------------------------------------
// Multi-Level DWT:
// For additional levels, use the LL subband as the new input and process similarly.
// -------------------------------------------------

// -------------------------------------------------
// Floating-Point Inverse DWT for Reconstruction:
// -------------------------------------------------
// Convert fixed-point subbands back to floating point using their LSB weights.
// Perform the inverse DWT (using standard floating-point operations) to reconstruct the image.

// -------------------------------------------------
// PSNR Calculation:
MSE = mean((original_image - reconstructed_image)^2);
PSNR = 10 * log10((255.0 * 255.0) / MSE);

if (PSNR >= 40)
    // Design meets specification
else
    // Adjust word lengths/scaling factors accordingly
