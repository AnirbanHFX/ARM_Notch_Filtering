; Owner
; Anirban Majumder
; Git : https://github.com/AnirbanHFX
; Provided as is

	AREA main, CODE, READONLY
	EXPORT __main
	ENTRY
	
inp dcw 0x0158, 0x0147, 0x003b, 0xfcf3, 0xfe21, 0xff7a, 0xfd69, 0xfb9b, 0xfbc5, 0xfe9f, 0xff26, 0xfe6d, 0xfb60, 0xfbcd, 0xfc10, 0xfa5c, 0xf8b8, 0xfc4b, 0x00e3, 0x0315, 0x022a, 0x0400, 0x0572, 0x0400, 0x0209, 0xfafb, 0xff04, 0xfc75, 0xfa08, 0xf7d5, 0xf94f, 0xf9d5, 0xf9f7, 0xf7bc, 0xf992, 0xfd71, 0xfb36, 0xf968, 0xfa53, 0xfada, 0xf4af, 0xf760, 0xf829, 0xfc4b, 0xf90c, 0xf34e, 0xf3bc, 0xf4ea, 0xf981, 0xf5ff, 0xf5dd, 0xf7a3, 0xf779, 0xf5dd, 0xf581, 0xfef3, 0x0305, 0xfdbd, 0xf9f7, 0xf821, 0xf91d, 0xf6af, 0xf45b, 0xf91d, 0xfce2, 0xfc08, 0xfa86, 0xfdb4, 0xfe65, 0xfe54, 0x003b, 0xfe9f, 0xfe86, 0xfc86, 0xf8fb, 0xfdbd, 0xff50, 0xfe3b, 0xfb82, 0xfc21, 0xfd60, 0x006d, 0xff36, 0xfd69, 0x00fc, 0x025d, 0xfe19, 0xfba3, 0xfa6d, 0xfa32, 0xf189, 0xf32d, 0xf4ea, 0xf5ee, 0xf979, 0xf947, 0xf93e, 0xfdd6, 0xfe11, 0xfe00, 0xff7a, 0x01b5, 0x047e, 0x0211, 0x012e, 0x04a0, 0x0326, 0xff60, 0xfa97, 0xfeda, 0x0115, 0xff58, 0x02b9, 0x0019, 0x0369, 0xff8a, 0xff2e, 0x0169, 0x0350, 0x0305, 0xff36, 0x0097, 0x0465, 0x0171, 0xffde, 0xfe19, 0xfcda, 0x01ac, 0x03ce, 0xfd3f, 0x00b9, 0x0019, 0xff50, 0xfbe6, 0xffa4, 0x0348, 0x03ce, 0x0393, 0x058b, 0x0812, 0x04a8, 0xfe11, 0xfb0c, 0xff60, 0xff15, 0xfd58, 0xfd58, 0xfc00, 0xfab0, 0xf9c5, 0xf71c, 0xfa5c, 0xfb79, 0xf903, 0xf7a3, 0xfc7e, 0xfe9f, 0xfc9f, 0xfdcd, 0xfec9, 0x0076, 0xfec1, 0xff36, 0x0287, 0x0315, 0x02c1, 0x045d, 0x0326, 0x0572, 0x03a4, 0x0498, 0x071e, 0x0a01, 0x04ec, 0x0076, 0x0011, 0x0411, 0x033f, 0xfd04, 0xfd69, 0xfde7, 0x00f4, 0x00d2, 0x0019, 0xffbd, 0x0126, 0x0000, 0xff50, 0x0422, 0x031e, 0x01df, 0x01ac, 0x046e, 0x0794, 0x084c, 0x010d, 0xffac, 0x038b, 0x00ca, 0x01e7, 0xfe54, 0x0315, 0x028f, 0x0147, 0x00a0, 0xffde, 0x01e7, 0xfe65, 0xfd71, 0xff36, 0xfe9f, 0xff58, 0xfe3b, 0xff36, 0x0008, 0xfb58, 0xfd15, 0xfd71, 0x023b, 0x0032, 0x01f8, 0x02a8, 0x0561, 0x0444, 0xffde, 0xffc5, 0x006d, 0x022a, 0x0000, 0x02ca, 0x03ac, 0x038b, 0xfeb0, 0xff71, 0xfea8, 0xfe6d, 0xfb47, 0xfbf7, 0xff1d, 0x02f4, 0x0104, 0x00b9, 0x00eb, 0x0772, 0x03e7, 0x02a8, 0x0032, 0x00d2, 0x00e3, 0xfd1d, 0xfd58, 0xfe4b, 0x0043, 0xfd60, 0xfea8, 0xfc00, 0x00da, 0xfa4b, 0xfe43, 0xff2e, 0x0076, 0xfc21, 0xf832, 0xfc97, 0xfdac, 0xfd25, 0xfd69, 0xfc64, 0xff15, 0xfda3, 0xfa8e, 0xfd47, 0xfad1, 0xfbef, 0xf9de, 0xfa5c, 0xfba3, 0xfbcd, 0xf9f7, 0xf8fb, 0xf8ea, 0xf9bc, 0xf725, 0xf420, 0xf9ef, 0xf90c, 0xf714, 0xf86c, 0xf6c0, 0xf83a, 0xf87d, 0xf618, 0xf8b8, 0xf925, 0xf818, 0xf3e6, 0xf8e2, 0xfc21, 0xfd9b, 0xf9e6, 0xf832, 0xfd25, 0xfd47, 0xf9cd, 0xf829, 0xf7cd, 0xf607, 0xf535, 0xf5e6, 0xf67d, 0xf864, 0xf83a, 0xf78a, 0xfa08, 0xff7a, 0xfc21, 0xf9cd, 0xfe11, 0x00e3, 0xfe43, 0xf8b8, 0xf8e2, 0xfd82, 0xff8a, 0xfdd6, 0xff69, 0x017a, 0x012e, 0x01df, 0x02d2, 0x02a0, 0xfec1, 0xfbe6, 0xf853, 0xfcc1, 0x0019, 0xfeda, 0xfd04, 0xfed2, 0x010d, 0xfe9f, 0xfc08, 0xfe11, 0xfd82, 0xfd47, 0xfa7d, 0xfe86, 0x0065, 0x00e3, 0xfeb0, 0xfe97, 0x0193, 0xff1d, 0xfeb9, 0xfe8f, 0x0126, 0x0211, 0xfd47, 0xfaea, 0xfca8, 0xfdef, 0xfc19, 0xf9bc, 0xfc10, 0x0022, 0xfbf7, 0xf8e2, 0xfc4b, 0xfd58, 0xfc19, 0xf79a, 0xfae2, 0xfc32, 0xffb4, 0xfe86, 0xfdd6, 0x0233, 0xffbd, 0xfd3f, 0xfd0c, 0x02a0, 0x01ac, 0x02a0, 0x02eb, 0x0444, 0x03ce, 0x0444, 0x02c1, 0x053f, 0x0750, 0x02f4, 0x0126, 0x0422, 0x0561, 0x0569, 0x0315, 0x041a, 0x0794, 0x06c2, 0x0422, 0x0158, 0x0759, 0x0716, 0x0305, 0xff04, 0x0065, 0x0097, 0x00a8, 0xfe75, 0xff82, 0xff71, 0xfc32, 0xfab0, 0xfe11, 0x02a8, 0x0287, 0x0008, 0x0382, 0x0411, 0x0361, 0xffce, 0xffbd, 0x0011, 0xfcda, 0xf9de, 0x0000, 0xffde, 0x0326, 0x02eb, 0x0348, 0x07f0, 0x0727, 0x038b, 0x0233, 0x0526, 0x09b5, 0x0750, 0x0444, 0x0801, 0x0a77, 0x0665, 0x0644, 0x086e, 0x0855, 0x06ca, 0x033f, 0x039c, 0x05f0, 0x05bd, 0x0372, 0x062b, 0x08b1, 0x0740, 0x0305, 0x0476, 0x058b, 0x084c, 0x05bd, 0x03b5, 0x06a0, 0x0887, 0x07be, 0x0348, 0x03a4, 0x04c2, 0xfc8e, 0x057a, 0x05f0, 0x078b, 0x058b, 0x01df, 0x02d2, 0x0569, 0x041a, 0x0043, 0x00da, 0x02b9, 0x04a8, 0x0369, 0x018b, 0x0382, 0x0644, 0x052f, 0x02d2, 0x0254, 0x0337, 0x033f, 0x0297, 0x032e, 0x028f, 0x0022, 0xff50, 0x02e3, 0x04db, 0x0801, 0x0561, 0x053f, 0x0905, 0x0372, 0x02ca
h	dcw 0x0000, 0xffe6, 0x0000, 0x0000, 0x0016, 0xfff8, 0x0000, 0x0000, 0x0000, 0x0000, 0xfff8, 0xffe5, 0x0000, 0x0017, 0x0016, 0x0000, 0x00e5, 0x0000, 0x0016, 0x0017, 0x0000, 0xffe5, 0xfff8, 0x0000, 0x0000, 0x0000, 0x0000, 0xfff8, 0x0016, 0x0000, 0x0000, 0xffe6, 0x0000
f equ 8			; Length of fractional part of fixed point
	
l equ 1000 		; Input length * 2 for 2 byte alignment  (500*2)
hl equ 66 		; Filter length * 2 for 2 byte alignment (33*2)
ol equ 1064 	; Output length * 2 for 2 byte alignment (532*2)

__mul		; Multiply two 16 bit fixed point numbers in r11, r12 and return in r10
			; Final result should be representable in A(8,8) fixed point
	mul r10, r11, r12
	asr r10, #f
	cmp r10, #0
	bmi __neg
	movt r10, #0x0000	; Extend upper 16 bits if number is positive
	bx lr
__neg
	movt r10, #0xFFFF	; Extend upper 16 bits if number is negative
	bx lr

__main
	
	ldr r0, =inp		; Address of input
	ldr r1, =h			; Address of filter coefficients
	mov r2, #0x0000		
	movt r2, #0x2000	; Address of output
	
	mov r3, #0		; i = 0
	
__outer_loop		; for(i=0; i<ol; i++)

	cmp r3, #ol
	bge __exit_outer	; exit if i >= ol
	mov r4, #0			; j = 0
	mov r9, #0			; accumulator = 0
	
__inner_loop		; for(j=0; j<hl; j++)

	cmp r4, #hl			; exit if j >= hl
	bge __exit_inner
	
	subs r5, r3, r4		; r5 = i - j
	
	bmi __skip			; Skip if i-j < 0
	cmp r5, #l			; Skip if i-j >= len(input)
	bge __skip
	
	add r6, r5, r0		; r6 = &x[i-j]
	ldrsh r11, [r6]		; r11 = x[i-j]
	add r7, r4, r1		; r7 = &h[j]
	ldrsh r12, [r7]		; r12 = h[j]
	bl __mul			; r10 = h[j]*x[i-j]
	add r9, r9, r10		; accumulator += h[j]*x[i-j]
	
__skip
	
	add r4, r4, #2
	b __inner_loop
__exit_inner

	add r8, r3, r2		; r8 = &y[i]
	strh r9, [r8]		; y[i] = accumulator

	add r3, r3, #2
	b __outer_loop
__exit_outer
	
	nop
	end