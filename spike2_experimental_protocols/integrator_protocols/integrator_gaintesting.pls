;-----------------------------------------------------------------------------
; Initial Parameters
;-----------------------------------------------------------------------------
            ; Set rate to 1 ms/step, 1 scaling, and 0 offset            
            SET    1.000,1,0

            ; Flags (1=running, 0=stopped)
            VAR    V1,BlockFlg=0   ;General block flag
            VAR    V2,Counter1     ;Main counter
            VAR    V3,Counter2     ;Sub-counter 1
            VAR    V4,Counter3     ;Sub-counter 2
            VAR    V5,Counter4     ;Sub-counter 3
            VAR    V6,NumGap=0
            VAR    V7,NumTest=0
            VAR    V8,NumTrain=0

            ; Drum parameters
            VAR    V9,DrumDac=0
            VAR    V10,DrumAmp=0
            VAR    V11,DrumOff=0

            ; Chair parameters
            VAR    V12,ChairAmp=0
            VAR    V13,ChairOff=0

            ; Sinusoidal stimulus parameters
            VAR    V14,SinFreq=0
            VAR    V15,SinPhase=0
            VAR    V16,SinAngle=0


;-----------------------------------------------------------------------------
; IDLELOOP: Sequencer loop when all signals are OFF
;-----------------------------------------------------------------------------
IDLELOOP: 'I DAC   0,0             ;Idling loop        >IDLING
            DAC    1,0             ;Repeatedly sets offset of drum and chair >"
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
            WAVEST S
            RATE   0,0
            RATE   1,0
            JUMP   IDLELOOP

;-----------------------------------------------------------------------------
; DSINE: Sinusoidal drum stimulus cycles
;-----------------------------------------------------------------------------
DSINE:  'D  SZ     0,DrumAmp       ;Start drum sine
            OFFSET 0,DrumOff       ;Set cosine offset
            PHASE  0,SinPhase      ;Set cosine phase
            ANGLE  0,SinAngle      ;Set cosine angle
            RATE   0,SinFreq       ;Set cosine frequency
DSINE1:     WAITC  0,DSINE1
            OFFSET 0,DrumOff
            JUMP   DSINE1          ;Set cosine loop    > Sine running

DSINEOFF: 'd CLRC  0               ;Stop drum sine at phase 0
DSINE2:     OFFSET 0,DrumOff       ;Adjust cosine offset
            WAITC  0,DSINE2        ;Wait for end of last cycle
            RATE   0,0,IDLELOOP    ;Stop chair sine then idle

;-----------------------------------------------------------------------------
; CSINE: Sinusoidal chair stimulus cycles
;-----------------------------------------------------------------------------
CSINE:  'C  SZ     1,ChairAmp      ;Start chair sine
            OFFSET 1,ChairOff      ;Set cosine offset
            PHASE  1,SinPhase      ;Set cosine phase
            ANGLE  1,SinAngle      ;Set the cosine angle
            RATE   1,SinFreq       ;Set cosine frequency
CSINE1:     WAITC  1,CSINE1
            OFFSET 1,ChairOff
            JUMP   CSINE1          ;Set cosine loop    > Sine running
CSINEOFF: 'c CLRC  1               ;Stop chair sine at phase 0
CSINE2:     OFFSET 1,ChairOff      ;Adjust cosine offset
            WAITC  1,CSINE2        ;Wait for end of last cycle
            RATE   1,0,IDLELOOP    ;Stop chair sine then idle

;-----------------------------------------------------------------------------
; SINE: Single Sinusoidal drum stimulus cycle
;-----------------------------------------------------------------------------
SINE:   'S  MOVI   BlockFlg,1      ;Start drum sine
            MOVI   Counter1,10
            SZ     1,ChairAmp      ;Set cosine amplitude
            OFFSET 1,ChairOff      ;Set cosine offset
            PHASE  1,SinPhase      ;Set cosine phase
            ANGLE  1,SinAngle      ;Set the cosine angle
            RATE   1,SinFreq       ;Set cosine frequency
SINE1:      WAITC  1,SINE1
            DBNZ   Counter1,SINE1  ;Run steps until counter hits zero >"
SINEOFF: 's CLRC   1               ;Stop chair sine at phase 0
SINE2:      WAITC  1,SINE2         ;Wait for end of last cycle
            RATE   1,0             ;Stop chair sine
            MOVI   BlockFlg,0      ;Set block as inactive
            JUMP   IDLELOOP

;-----------------------------------------------------------------------------
; OCULOMOTOR: Oculomotor Integrator stimulus
;-----------------------------------------------------------------------------
OCULOON: 'O DAC    0,DrumDac       ;Start oculumotor integrator
            JUMP   OCULOON         ;Loop               > Oculo running

OCULOOFF: 'o DAC   0,0             ;Stop oculomotor integrator
            DAC    1,0             ;Reset chair offset
            JUMP   IDLELOOP

;-----------------------------------------------------------------------------
; GAP: No stimuli for Gap1Dur+Gap3Dur ms and FlashDur ms light pulse halfway
;-----------------------------------------------------------------------------
GAP:    'G  MOVI   BlockFlg,1      ;Start Gap block    >GAP BLOCK
            DIGOUT [.......0]      ;Reset to initial state >"
            MOV    Counter1,NumGap ;Set number of step periods to run >"
GAP1:       DBNZ   Counter1,GAP1   ;Run steps until counter hits zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP

;-----------------------------------------------------------------------------
; TRAIN: x2/x1/x0 training block
;-----------------------------------------------------------------------------
TRAIN:  'T  MOVI   BlockFlg,1      ;Start training block >TRAINING
            DIGOUT [.......0]      ;Ensure light is off >"
            MOV    Counter3,NumTrain ;Set number of step periods to run >"
            DIGOUT [.......1]      ;Turn light on      >"
TRAIN1:     DAC    0,DrumDac
            DBNZ   Counter3,TRAIN1 ;Run steps until counter hits zero >"
            DIGOUT [.......0]      ;Turn light back off >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TEST: pre/post test VORD (chair only) block
;-----------------------------------------------------------------------------
TEST:   'P  MOVI   BlockFlg,1      ;Start Test block   >TESTING
            DIGOUT [.......0]      ;Ensure light is off >"
            MOV    Counter2,NumTest ;Set number of step periods to run >"
TEST1:      DBNZ   Counter2,TEST1  ;Run steps until counter hits zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP

;-----------------------------------------------------------------------------
; WAVE: Custom waveform stimuli (if provided)
;-----------------------------------------------------------------------------
;WAVE:   'W  WAVEGO w,TW
;            WAVEST T