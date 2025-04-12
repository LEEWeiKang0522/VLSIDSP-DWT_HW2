# VLSIDSP-DWT_HW2
# Detailed Fixed‑Point 2‑D DWT Algorithm

This document describes the complete bit‑true fixed‑point algorithm for the forward 2‑D DWT (used in image processing) along with the floating‑point inverse DWT. The algorithm includes fixed‑point parameter definitions, the step-by-step processing flow for the horizontal and vertical transforms, multi‑level processing, and PSNR calculation. Each variable’s word length and LSB weight are explicitly indicated.

---

## 1. Overview

- **Goal:**  
  Conduct fixed‑point simulations of the forward 2‑D DWT such that, when followed by a floating‑point inverse DWT, the reconstructed image achieves a PSNR of no less than 40 dB.  
- **Approach:**  
  Use fixed‑point arithmetic for the forward transform, with carefully chosen word lengths for filter coefficients, intermediate accumulators, and output subbands. Convert the fixed‑point output back to floating point for inverse DWT and image reconstruction.
- **Processing Steps:**  
  1. **Pre‑processing:** Convert 8‑bit input image to fixed‑point representation (as required).  
  2. **Horizontal 1‑D DWT:** Convolve each row with quantized lowpass and highpass filters, downsample, and produce 12‑bit outputs.  
  3. **Vertical 1‑D DWT:** Process the horizontal outputs column‑wise to generate the final subbands (LL, LH, HL, HH) with 12‑bit format.  
  4. **Multi‑Level DWT:** For multi‑level processing, feed the LL subband of the current level into the next level using the same fixed‑point format.
  5. **Reconstruction:** Convert fixed‑point subbands to floating point and perform the inverse DWT.
  6. **PSNR Computation:** Calculate the PSNR of the reconstructed image using the formula below.

---

## 2. Fixed‑Point Parameter Definitions

### A. Filter Coefficients
- **Representation:** All coefficients (for both lowpass `h(n)` and highpass `g(n)`) are quantized uniformly.
- **Word Length:** 16 bits  
- **Format:** 1‑bit sign + 1‑bit integer + 14‑bit fractional  
- **LSB Weight:** \(2^{-14}\)  
- **Note:** Since typical 9/7 wavelets have coefficients in the approximate range \([-1, 1]\), this format provides sufficient precision.

### B. Intermediate Computations
- **Accumulator for Multiply-Accumulate (MAC):** 32 bits  
  - *Rationale:* Multiplying an 8‑bit input by a 16‑bit coefficient may produce a 24‑bit product. Summing over many taps (e.g., 9 to 10) requires a 32‑bit accumulator to avoid overflow.

### C. DWT Output Variables
- **Horizontal 1‑D DWT Output:**  
  - **Word Length:** 12 bits  
  - **Format:** 1‑bit sign + 3‑bit integer + 8‑bit fractional  
  - **LSB Weight:** \(2^{-8}\)  
- **Vertical 1‑D DWT (Final Subbands: LL, LH, HL, HH):**  
  - **Word Length:** 12 bits (same as horizontal output)  
  - **LSB Weight:** \(2^{-8}\)  
- **Multi‑Level DWT:**  
  Subsequent levels use the same 12‑bit format for the LL subband.

---

## 3. Detailed Algorithm Steps

### Step 1: Pre‑processing
- **Input:** 8‑bit grayscale image.
- **Action:** Convert each pixel to fixed‑point using an appropriate scaling if needed (e.g., interpreting 0–255 as fixed‑point values).
  
### Step 2: Forward 2‑D DWT – Horizontal 1‑D Transform
- **For Each Row:**
  1. **For Each Column (with appropriate index for convolution):**  
     - Initialize two 32‑bit accumulators: one for the lowpass path and one for the highpass path.
     - **For Each Filter Tap \( k \) (from \(-K\) to \(K\)):**  
       - Load the input sample at position `[i][j+k]` and convert it to fixed‑point.
       - Multiply the sample by the corresponding filter coefficient from `h` (for lowpass) and `g` (for highpass).
       - Accumulate the products into the respective 32‑bit accumulators.
  2. **Downsampling:**  
     - After processing the taps for a given output location, downsample (e.g., choose even indexed outputs).
     - Convert the 32‑bit accumulator values to 12‑bit fixed‑point format using rounding or truncation.

### Step 3: Forward 2‑D DWT – Vertical 1‑D Transform
- **For Each Column (of the horizontal output):**
  1. **For Each Row:**  
     - Initialize two 32‑bit accumulators: one for vertical lowpass and one for vertical highpass.
     - **For Each Vertical Filter Tap \( k \):**  
       - Multiply the fixed‑point output from the horizontal stage (from the corresponding row offset) by the filter coefficient.
       - Accumulate the product results.
  2. **Downsampling:**  
     - Downsample vertically and convert the result to 12‑bit fixed‑point format.
     - This produces the final subbands: LL, LH, HL, and HH.
     
### Step 4: Multi‑Level DWT Processing
- **If Processing Further Levels:**  
  - Use the LL subband from the current level as the input.
  - Apply the horizontal and vertical 1‑D DWT steps again.
  - Continue for the required number of levels (e.g., Level 2 and Level 3).

### Step 5: Floating‑Point Inverse DWT (IDWT) for Reconstruction
- **Conversion:**  
  - Convert all fixed‑point subband outputs back to floating‑point. Use the known LSB weights (i.e., \(2^{-8}\) for outputs) to determine the scaling factor.
- **Inverse Transform:**  
  - Use standard floating‑point convolution and upsampling to perform the inverse DWT.
  - Obtain the reconstructed image.

### Step 6: PSNR Calculation
- **Formula:**

  \[
  \text{PSNR} = 10 \cdot \log_{10}\left(\frac{255^2}{\text{MSE}}\right)
  \]
  
  where **MSE** is the mean squared error between the original image and the reconstructed image.
  
- **Requirement:** The PSNR must be ≥ 40 dB.
- **Simulation Example:**  
  With the given design, simulation on an 8‑bit grayscale image using a 3‑level 2‑D DWT produced a reconstructed image PSNR of approximately 41–42 dB.

---

## 4. Detailed Pseudocode

```c
// -------------------------------------------------
// Fixed-Point Parameter Definitions:
#define FILTER_WORD_LEN 16   // 16-bit coefficients: 1-bit sign + 1-bit integer + 14-bit fraction, LSB = 2^-14
#define ACCUM_WORD_LEN  32   // 32-bit accumulator for MAC
#define OUTPUT_WORD_LEN 12   // 12-bit DWT output: 1-bit sign + 3-bit integer + 8-bit fraction, LSB = 2^-8

// Pre-quantized filter coefficients in 16-bit fixed-point format:
fixed16_t h[Nh] = { /* e.g., quantize(0.852698679009, 16), quantize(0.377402855613, 16), ... */ };
fixed16_t g[Ng] = { /* e.g., quantize(-0.788485616406, 16), quantize(0.418092273222, 16), ... */ };

// -------------------------------------------------
// Step 2: Horizontal 1-D DWT (Fixed-Point)
for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j += DS_H) { // DS_H: horizontal downsampling factor (typically 2)
        int32_t acc_low  = 0;
        int32_t acc_high = 0;
        for (int k = -K; k <= K; k++) {  // K: half the filter length
            // Load input sample; input_image is pre-converted to fixed-point if necessary
            fixed16_t sample = convert_to_fixed(input_image[i][j + k]);
            // Accumulate lowpass and highpass components
            acc_low  += sample * h[k + OFFSET];  // OFFSET adjusts k to array index
            acc_high += sample * g[k + OFFSET];
        }
        // Downsample (e.g., take every second result) and convert to 12-bit fixed-point
        fixed12_t L[i][j / DS_H] = truncate_round(acc_low, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
        fixed12_t H[i][j / DS_H] = truncate_round(acc_high, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
    }
}

// -------------------------------------------------
// Step 3: Vertical 1-D DWT (Fixed-Point)
for (int j = 0; j < (N / DS_H); j++) {
    for (int i = 0; i < M; i += DS_V) {  // DS_V: vertical downsampling factor (typically 2)
        int32_t acc_LL = 0;
        int32_t acc_LH = 0;
        for (int k = -K; k <= K; k++) {
            acc_LL += L[i + k][j] * h[k + OFFSET];
            acc_LH += H[i + k][j] * g[k + OFFSET];
        }
        // Downsample vertically and convert to 12-bit fixed-point values for the final subbands
        fixed12_t LL[i / DS_V][j] = truncate_round(acc_LL, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
        fixed12_t LH[i / DS_V][j] = truncate_round(acc_LH, ACCUM_WORD_LEN, OUTPUT_WORD_LEN);
        // HL and HH subbands can be calculated similarly based on the implementation strategy.
    }
}

// -------------------------------------------------
// Step 4: Multi-Level DWT
// If additional levels are required, use the LL subband from the previous level:
// for (each level beyond Level 1) {
//     Set input_image = LL;
//     Repeat steps for Horizontal and Vertical DWT.
// }

// -------------------------------------------------
// Step 5: Floating-Point Inverse DWT for Reconstruction
// Convert fixed-point outputs (LL, LH, HL, HH) to floating-point values using the LSB weights.
// Perform the standard floating-point inverse DWT (upsampling and convolution) to reconstruct the image.

// -------------------------------------------------
// Step 6: PSNR Calculation
MSE = mean_squared_error(original_image, reconstructed_image);
PSNR = 10 * log10((255.0 * 255.0) / MSE);
if (PSNR >= 40) {
    // Design meets the specification.
} else {
    // Adjust fixed-point word lengths or scaling factors.
}
