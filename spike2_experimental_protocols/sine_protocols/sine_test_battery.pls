;-----------------------------------------------------------------------------
; SINE - TEST BATTERY PROTOCOL OUTPUT SEQUENCER SCRIPT (SPIKE2 VERSION 10)
;-----------------------------------------------------------------------------
; Initial Parameters
;-----------------------------------------------------------------------------
            ; Set rate to 1 ms/sequencer step, 1 scaling, and 0 offset
            SET    1.000,1,0

            ; --- Core Sequencer Variables ---
            ; These variables should be the same across ALL protocols
            ; Note V1-V256 is available, but DO NOT use V56-V64, V255, V256

            ; Block Flag V1 (tells Spike2 if a block is actively running)
            VAR    V1,BlockFlg=0   ;1 = running, 0 = not running

            ; Drum/Target parameters (V2 - V9)
            VAR    V2,DrumOff=0    ;Drum offset
            VAR    V3,DrumAmp=0    ;Drum amplitude
            VAR    V4,DrumFrq=0    ;Drum frequency
            VAR    V5,DrumPhs=0    ;Drum phase
            VAR    V6,DrumAng=0    ;Drum angle
            VAR    V7,NdrumP=0     ;Number of drum TEST block cycles
            VAR    V8,NdrumT=0     ;Number of drum TRAIN block cycles
            VAR    V9,DrumCtr=0    ;Cycle counter for drum

            ; Chair/Head parameters (V10 - V17)
            VAR    V10,ChairOff=0  ;Chair offset
            VAR    V11,ChairAmp=0  ;Chair amplitude
            VAR    V12,ChairFrq=0  ;Chair frequency
            VAR    V13,ChairPhs=0  ;Chair phase
            VAR    V14,ChairAng=0  ;Chair angle
            VAR    V15,NchairP=0   ;Number of chair TEST block cycles
            VAR    V16,NchairT=0   ;Number of chair TRAIN block cycles
            VAR    V17,ChairCtr=0  ;Cycle counter for chair

            ; Gap block timing parameters (V18 - V21)
            VAR    V18,Gap1Dur=478 ;(Gap1Dur*rate/5)-2
            VAR    V19,FlashDur=39 ;(FlashDur*rate/5)-1
            VAR    V20,Gap3Dur=479 ;(Gap3Dur*rate/5)-1
            VAR    V21,GapCtr=0    ;Counter for gap block

            ; Variables for storing dynamic drum/chair amp values
            VAR    V100,DrumTmp=0  ;Variable drum amplitude
            VAR    V101,ChairTmp=0 ;Variable drum amplitude


;-----------------------------------------------------------------------------
; IDLELOOP: Sequencer loop when all signals are OFF
;-----------------------------------------------------------------------------
IDLELOOP: 'I DAC   0,DrumOff       ;Idling loop        >IDLING
            DAC    1,ChairOff      ;Repeatedly sets offset of drum and chair >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TTL1OFF: Turns off TTL channel 1 (digital output bit 8)
;-----------------------------------------------------------------------------
TTL1OFF: 'l DIGOUT [.......0]      ;Turn TTL 1 off     >=
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TTL1ON: Turns on TTL channel 1 (digital output bit 8)
;-----------------------------------------------------------------------------
TTL1ON: 'L  DIGOUT [.......1]      ;Turn TTL 1 on      >=
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; RESET: Resets sequencer to initial state
;-----------------------------------------------------------------------------
RESET:  'R  DIGOUT [.......0]      ;Reset to initial state >=
            RATE   0,0             ;Stop sine on drum  >=
            RATE   1,0             ;Stop sine on chair >=
            JUMP   IDLELOOP        ;Return to idle loop


;-----------------------------------------------------------------------------
; DRUM SINE: Turn sinusoidal drum only stimulus on/off
;-----------------------------------------------------------------------------
DSINEON: 'D MOV    DrumTmp,DrumAmp ;Set drum amplitude >"
            MOVI   ChairTmp,0      ;Set chair amplitude >"
            SZ     0,DrumAmp       ;Start drum cosine >"
            OFFSET 0,DrumOff       ;Set cosine offset >"
            PHASE  0,DrumPhs       ;Set cosine relative phase >"
            ANGLE  0,0             ;Set cosine angle >"
            RATE   0,DrumFrq       ;Set cosine frequency >"
DSINE1:     WAITC  0,DSINE1        ;Wait for 0 drum phase >"
            OFFSET 0,DrumOff       ;                   >"
            JUMP   DSINE1          ;Set cosine loop    > Drum running

DSINEOFF: 'd CLRC  0               ;Stop drum sine at phase 0 >"
DSINE2:     OFFSET 0,DrumOff       ;Adjust cosine offset >"
            WAITC  0,DSINE2        ;Wait for end of cycle >"
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            RATE   0,0,IDLELOOP    ;Stop drum cosine then idle


;-----------------------------------------------------------------------------
; CHAIR SINE: Turn sinusoidal chair only stimulus on/off
;-----------------------------------------------------------------------------
CSINEON: 'C MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 1,ChairOff      ;Set cosine offset >"
            PHASE  1,ChairPhs      ;Set cosine relative phase >"
            ANGLE  1,0             ;Set cosine angle >"
            RATE   1,ChairFrq      ;Set cosine frequency >"
CSINE1:     WAITC  1,CSINE1        ;Wait for 0 chair phase >"
            OFFSET 1,ChairOff      ;                   >"
            JUMP   CSINE1          ;Set cosine loop    > Chair running

CSINEOFF: 'c CLRC  1               ;Stop chair sine at phase 0 >"
CSINE2:     OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,CSINE2        ;Wait for end of cycle >"
            MOVI   ChairTmp,0      ;Set chair amplitude >"
            RATE   1,0,IDLELOOP    ;Stop chair cosine


;-----------------------------------------------------------------------------
; SINE: Turn sinusoidal drum and chair stimulus on/off
;-----------------------------------------------------------------------------
SINEON: 'S  MOV    DrumTmp,DrumAmp ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            SZ     0,DrumAmp       ;Start sine stimulus >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 0,DrumOff       ;Set drum cosine offset >"
            OFFSET 1,ChairOff      ;Set chair cosine offset >"
            PHASE  0,DrumPhs       ;Set drum cosine relative phase >"
            PHASE  1,ChairPhs      ;Set chair cosine relative phase >"
            ANGLE  0,0             ;Set drum cosine angle >"
            ANGLE  1,0             ;Set chair cosine angle >"
            RATE   0,DrumFrq       ;Set drum cosine frequency >"
            RATE   1,ChairFrq      ;Set chair cosine frequency >"
SINE1:      OFFSET 0,DrumOff       ;Adjust drum cosine offset >"
            OFFSET 1,ChairOff      ;Adjust chair cosine offset >"
            WAITC  1,SINE1         ;Wait for 0 chair phase >"
            JUMP   SINE1           ;Set cosine loop    > Sine running

SINEOFF: 's CLRC   1               ;Stop sines at phase 0 of chair
SINE2:      OFFSET 0,DrumOff       ;Adjust drum cosine offset >"
            OFFSET 1,ChairOff      ;Adjust chair cosine offset >"
            WAITC  1,SINE2         ;Wait for end of cycle >"
            RATE   0,0             ;Stop drum cosine   >"
            RATE   1,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOVI   ChairTmp,0      ;Set chair amplitude >"
            JUMP   IDLELOOP        ;Return to idle loop


;-----------------------------------------------------------------------------
; GAP: No stimuli for Gap1Dur+Gap3Dur ms and FlashDur ms light pulse halfway
;-----------------------------------------------------------------------------
GAP:    'G  MOVI   BlockFlg,1      ;Start Gap block    >GAP BLOCK
            DIGOUT [.......0]      ;Reset to initial state >"
            RATE   0,0             ;Stop cosine on drum >"
            RATE   1,0             ;Stop cosine on chair >"
            DAC    0,DrumOff       ;Stop the drum      >"
            DAC    1,ChairOff      ;Stop the chair     >"
            MOV    GapCtr,Gap1Dur  ;Set duration of first half >"
GAP1:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   GapCtr,GAP1     ;Repeat until counter hits zero >"
            BEQ    FlashDur,0,GSKIP ;Skip to GAP2 if FlashDur is zero >"
            DIGOUT [.......1]      ;Turn light on      >"
            MOV    GapCtr,FlashDur ;Set duration of light pulse >"
GAP2:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   GapCtr,GAP2     ;Repeat until counter hits zero >"
GSKIP:      DIGOUT [.......0]      ;Turn light off     >"
            MOV    GapCtr,Gap3Dur  ;Set duration of second half >"
GAP3:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   GapCtr,GAP3     ;Repeat until counter hits zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; VORx0 Block (Light On, Drum & Chair In-Phase)
;-----------------------------------------------------------------------------
VOR0ON: 'V  MOVI   BlockFlg,1      ;Start VORx0 block  >VORx0
            NEG    DrumTmp,DrumAmp ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            MOV    ChairCtr,NchairP ;Set number of cycles to run >"
            DIGOUT [.......1]      ;Turn on light      >"
            SZ     0,DrumTmp       ;Start drum cosine  >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 0,DrumOff       ;Set drum cosine offset >"
            OFFSET 1,ChairOff      ;Set chair cosine offset >"
            PHASE  0,DrumPhs       ;Set drum cosine relative phase >"
            PHASE  1,ChairPhs      ;Set chair cosine relative phase >"
            ANGLE  0,0             ;Set drum cosine angle >"
            ANGLE  1,0             ;Set chair cosine angle >"
            RATE   0,DrumFrq       ;Set drum cosine frequency >"
            RATE   1,ChairFrq      ;Set chair cosine frequency >"
VOR01:      OFFSET 0,DrumOff       ;Adjust drum cosine offset >"
            OFFSET 1,ChairOff      ;Adjust chair cosine offset >"
            WAITC  1,VOR01         ;Wait for 0 chair phase >"
            DBNZ   ChairCtr,VOR01  ;Run cycles until counter hits zero >"

VOR0OFF: 'v CLRC   1               ;Stop chair sine at 0 phase >"
VOR02:      OFFSET 0,DrumOff       ;Adjust drum cosine offset >"
            OFFSET 1,ChairOff      ;Adjust chair cosine offset >"
            WAITC  1,VOR02         ;Wait for end of cycle >"
            RATE   0,0             ;Stop drum cosine   >"
            RATE   1,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude to zero >"
            MOVI   ChairTmp,0      ;Set chair amplitude to zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   TTL1OFF


;-----------------------------------------------------------------------------
; VORx1 Block (Light on, Chair only)
;-----------------------------------------------------------------------------
VOR1ON: 'W  MOVI   BlockFlg,1      ;Start VORx1 block  >VORx1
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            MOV    ChairCtr,NchairP ;Set number of cycles to run >"
            DIGOUT [.......1]      ;Turn on light      >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 1,ChairOff      ;Set cosine offset  >"
            PHASE  1,ChairPhs      ;Set cosine relative phase >"
            ANGLE  1,0             ;Set cosine angle   >"
            RATE   1,ChairFrq      ;Set cosine frequency >"
VOR11:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,VOR11         ;Wait for 0 phase   >"
            DBNZ   ChairCtr,VOR11  ;Run cycles until counter hits zero >"

VOR1OFF: 'w CLRC   1               ;Stop chair sine at 0 phase >"
VOR12:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,VOR12         ;Wait for end of cycle >"
            RATE   1,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude to zero >"
            MOVI   ChairTmp,0      ;Set chair amplitude to zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   TTL1OFF


;-----------------------------------------------------------------------------
; VORx2 Block (Light On, Chair & Drum Out-of-Phase)
;-----------------------------------------------------------------------------
VOR2ON: 'X  MOVI   BlockFlg,1      ;Start VORx2 block  >VORx2
            MOV    DrumTmp,DrumAmp ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            MOV    ChairCtr,NchairP ;Set number of cycles to run >"
            DIGOUT [.......1]      ;Turn on light      >"
            SZ     0,DrumTmp       ;Start drum cosine  >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 0,DrumOff       ;Set drum cosine offset >"
            OFFSET 1,ChairOff      ;Set chair cosine offset >"
            PHASE  0,DrumPhs       ;Set drum cosine relative phase >"
            PHASE  1,ChairPhs      ;Set chair cosine relative phase >"
            ANGLE  0,0             ;Set drum cosine angle >"
            ANGLE  1,0             ;Set chair cosine angle >"
            RATE   0,DrumFrq       ;Set drum cosine frequency >"
            RATE   1,ChairFrq      ;Set chair cosine frequency >"
VOR21:      OFFSET 0,DrumOff       ;Adjust drum cosine offset >"
            OFFSET 1,ChairOff      ;Adjust chair cosine offset >"
            WAITC  1,VOR21         ;Wait for 0 chair phase >"
            DBNZ   ChairCtr,VOR21  ;Run cycles until counter hits zero >"

VOR2OFF: 'x CLRC   1               ;Stop chair sine at 0 phase >"
VOR22:      OFFSET 0,DrumOff       ;Adjust drum cosine offset >"
            OFFSET 1,ChairOff      ;Adjust chair cosine offset >"
            WAITC  1,VOR22         ;Wait for end of cycle >"
            RATE   0,0             ;Stop drum cosine   >"
            RATE   1,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude to zero >"
            MOVI   ChairTmp,0      ;Set chair amplitude to zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   TTL1OFF


;-----------------------------------------------------------------------------
; VORD Block (Light off, Chair only)
;-----------------------------------------------------------------------------
VORDON: 'Y  MOVI   BlockFlg,1      ;Start VORD block   >VORD
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            MOV    ChairCtr,NchairP ;Set number of cycles to run >"
            DIGOUT [.......0]      ;Ensure light is off >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 1,ChairOff      ;Set cosine offset  >"
            PHASE  1,ChairPhs      ;Set cosine relative phase >"
            ANGLE  1,0             ;Set cosine angle   >"
            RATE   1,ChairFrq      ;Set cosine frequency >"
VORD1:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,VORD1         ;Wait for 0 phase   >"
            DBNZ   ChairCtr,VORD1  ;Run cycles until counter hits zero >"

VORDOFF: 'y CLRC   1               ;Stop chair sine at 0 phase >"
VORD2:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,VORD2         ;Wait for end of cycle >"
            RATE   1,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude to zero >"
            MOVI   ChairTmp,0      ;Set chair amplitude to zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; OKR Block (Light on, Drum only)
;-----------------------------------------------------------------------------
OKRON:  'Z  MOVI   BlockFlg,1      ;Start OKR block    >OKR
            MOV    DrumTmp,DrumAmp ;Set drum amplitude >"
            MOVI   ChairTmp,0      ;Set chair amplitude >"
            MOV    DrumCtr,NdrumP  ;Set number of cycles to run >"
            DIGOUT [.......1]      ;Turn light on      >"
            SZ     0,DrumAmp       ;Start chair cosine >"
            OFFSET 0,DrumOff       ;Set cosine offset  >"
            PHASE  0,DrumPhs       ;Set cosine relative phase >"
            ANGLE  0,0             ;Set cosine angle   >"
            RATE   0,DrumFrq       ;Set cosine frequency >"
OKR1:       OFFSET 0,DrumOff       ;Adjust cosine offset >"
            WAITC  0,OKR1          ;Wait for 0 phase   >"
            DBNZ   DrumCtr,OKR1    ;Run cycles until counter hits zero >"

OKROFF: 'z  CLRC   0               ;Stop chair sine at 0 phase >"
OKR2:       OFFSET 0,DrumOff       ;Adjust cosine offset >"
            WAITC  0,OKR2          ;Wait for end of cycle >"
            RATE   0,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude to zero >"
            MOVI   ChairTmp,0      ;Set chair amplitude to zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   TTL1OFF


;-----------------------------------------------------------------------------
; WAVE: Custom waveform stimuli (if provided)
;-----------------------------------------------------------------------------
;WAVE:   'W  WAVEGO w,TW
;            WAVEST T