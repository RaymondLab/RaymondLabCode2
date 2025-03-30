// matint.h
/*
    Copyright (C) Cambridge Electronic Design Limited 2014
    Author: James Thompson
    Web: www.ced.co.uk email: james@ced.co.uk, softhelp@ced.co.uk

    This file is part of CEDS64ML, a MATLAB interface to the SON64 library.

    CEDS64ML is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CEDS64ML is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CEDS64ML.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __MATINT_H__
#define __MATINT_H__

#define MAX_MATINT_FILES 100
#define MAX_MATINT_FILTERS 1000

#ifdef MATINT_EXPORTS
#define MATINT_API __declspec(dllexport)
#else
#define MATINT_API __declspec(dllimport)
#endif

#ifdef __cplusplus
extern "C" {
#endif
    /*Data Structures*/

    //! A structure for the marker data type
    /*!
    A marker consists of 64-bit time stamp and 4 8-bit chars that contain code infomation.
    */
    typedef struct S64Marker{
        long long m_Time;   //!< A 64-bit time stamp for the marker
        unsigned char m_Code1;  //!< The first marker code
        unsigned char m_Code2;  //!< The second marker code
        unsigned char m_Code3;  //!< The third marker code
        unsigned char m_Code4;  //!< The fourth marker code
    } S64Marker;

    /*Basic Functions*/

    //! Get the number of open files
    /*!
    Returns the number of SON files that this dll currently has open. (Unique pointers to) the files are
    held in a vector and at the moment we restrict the vector to length 8 so only 8 files can be open 
    simultaneously

    See also: S64Open(), S64IsOpen(), S64Create(), S64Close(), S64CloseAll()
    \return The number of open files.
    */
    MATINT_API int S64FileCount();
    //! Close all the currently open files
    /*!
    Closes all the currently open files

    See also: S64Open(), S64IsOpen(), S64Create(), S64FileCount(), S64CloseAll()
    \return The number of open files after attempting to close them all. If this number is non-zero then
    an error has occured.
    */
    MATINT_API int S64CloseAll();

    /*Opening/closing etc. functions*/

    //! Create a new, empty file on disk
    /*!
    Creates a new, empty file on disk. You will have to create any channels seperately with S64SetXChan()

    See also: S64Open(), S64IsOpen(), S64Close(), S64CloseAll(), S64SetXChan()
    \param FileName The path and name of the file. This follows the usual naming conventions of
    the system.
    \param nChans The initial number of channels. It is theoretically possible to increase the
    number of channels later, but there is not yet code to do this. The minimum
    number is 32. There is an upper limit of MAXCHANS (1000) but this is an
    arbitrary limit and could be increased. However, we only allow up to 128 header
    blocks, which imposes some limit on the absolute maximum number.
    \param nBig A flag to set the type of file. 0 = Very old 32-bit .smr format with a size limit of 2 GB, 
    1 = a 32-bit .smr file with as size limit of 1 TB, 2 = 64-bit .smrx file.
    This is only relevant when the file name does not contain a '.smr(x)' suffix
    \return If successful, a positive integer that will act as handle to this file in MATLAB.
    Otherwise a negative error code.
    */
    MATINT_API int S64Create(const char *FileName, const int nChans, const int nBig);

    //! Open an existing file
    /*!
    This command opens an existing data file for reading or reading and writing. If you open a
    file in read-only mode, you cannot use any commands that would change the data.

    See also: S64IsOpen(), S64Create(), S64Close(), S64CloseAll()
    \param FileName The path and name of the file. This follows the usual naming conventions of
    the system.
    \param nFlag 1=read only, 0=read/write, -1 readwrite, then try readonly
    \return If successful, a positive integer that will act as handle to this file in MATLAB.
    Otherwise a negative error code.
    */
    MATINT_API int S64Open(const char *FileName, const int nFlag);

    //! Check whether a file handle is currently assigned to an open file
    /*!
    This unction allows you to check if a file handle is currently in use, i.e. there is an open
    file with handle

    See also: S64Open(), S64Create(), S64Close(), S64CloseAll(), S64FileCount()
    \param nFid An integer file handle (1-8)
    An integer between 1 and 8 at the moment
    \return     1 if the file handle is in use, 0 if it is not in use, -1 if nFid is invalid
    */
    MATINT_API int S64IsOpen(const int nFid);

    //! Close an open file
    /*!
    Closes the file with file handle nFid

    See also: S64Open(), S64Create(), S64CloseAll(), S64FileCount()
    \param nFid An integer file handle (1-8)
    \return Returns 0 if the file was closed successfully or the was no file with that handle,
    otherwise returns a negative error message.
    */
    MATINT_API int S64Close(const int nFid);

    //! Empties a file of all data
    /*!
    This function should be called to clear all data from a file while leaving the channels intact
    It doesn't actually physically delete all the data from the disk, it just marks it as gone
    to allow it to be overwritten

    \param nFid An integer file handle (1-8)
    \return Returns 0 if the file was closed successfully or the was no file with that handle,
    otherwise returns a negative error message.
    */
    MATINT_API int S64Empty(const int nFid);

    /*Get/Set global properties for the file*/

    //! Get a file comment
    /*!
    Used to get a file comment. This function copies data from one char buffer to another. It is
    *CRITICAL* that the buffers are the right size or *MATLAB WILL CRASH*. What you have to do call this
    function with nGetSize set to a negative number. The function wont copy any data but will return
    the size of the required buffer. Then, in MATLAB, use this number to a create buffer of the correct
    size and call the function again with nGetSize set to a non-negative number
    See SetFileComment() for more information.
    \param nFid An integer file handle (1-8)
    \param nInd The index of the comment, 1 to 5
    \param Comment A string holding the comment to set.
    \param nGetSize A flag to determine the behaviour of the function. If ngetSize < 0, the function
    returns the length of the comment, but does not get the comment. If nGetSize >= 0 the function
    gets the comment
    \return If ngetSize < 0, the length of the comment of a negative error message
    If nGetSize >= 0, 0 or a negative code.
    */
    MATINT_API int S64GetFileComment(const int nFid, const int nInd, char *Comment, const int nGetSize);

    //! Set a file comment
    /*!
    There are up to NUMFILECOMMENTS (8) comments stored in the file header. There is
    no limit (save common sense and available space in the header) to the size of each
    comment, and comments can contain end of line characters - in fact anything that
    can be stored in a std::string of char. The 32-bit system limited comments to 79
    characters.
    \param nFid An integer file handle (1-8)
    \param nInd The index of the comment, 1 to 5
    \param Comment A string holding the comment to set.
    \return 0 or a negative code (BAD_PARAM if n is out of range).
    */
    MATINT_API int S64SetFileComment(const int nFid, const int nInd, const char *Comment);

    //! Find the lowest numbered channel that is not in use
    /*!
    \param nFid An integer file handle (1-8)
    \return The lowest channel number of a channel that is not used or error NO_CHANNEL.
    */
    MATINT_API int S64GetFreeChan(const int nFid);

    //! Get the maximum number of channels in the file
    /*!
    \param nFid An integer file handle (1-8)
    \return The maximum numberof channels that can be stored in the file.
    */
    MATINT_API int S64MaxChans(const int nFid);

    //! Get the seconds per clock tick
    /*!
    Everything in the file is quantified to the underlying clock tick. As all values in the
    file are stored, set and returned in ticks, you need to read this value to interpret
    times in seconds.
    \param nFid An integer file handle (1-8)
    \return The underlying clock tick period in seconds.
    */
    MATINT_API double S64GetTimeBase(const int nFid);

    //! Set the seconds per clock tick
    /*!
    \param nFid An integer file handle (1-8)
    \param dSecPerTick This sets the underlying clock tick period in seconds. Setting a
    value <= 0.0 has no effect.
    \return 0 or a negative code.
    */
    MATINT_API int S64SetTimeBase(const int nFid, const double dSecPerTick);

    //! Converts seconds to file ticks
    /*!
    \param nFid An integer file handle (1-8)
    \param dSec time in seconds

    \return dSec converted to time in ticks or a negative error code.
    */
    MATINT_API long long S64SecsToTicks(const int nFid, const double dSec);

    //! Converts file ticks to seconds
    /*!
    \param nFid An integer file handle (1-8)
    \param tSec time in seconds
    \return tSec converted to time in ticks or a negative error code.
    */
    MATINT_API double S64TicksToSecs(const int nFid, const long long tSec);

    //! Get the version of the data file
    /*!
    The returned value is the major file system version * 256 + the minor version. The
    major version number of the 32-bit filing system opened through this interface is 0.
    The first major version of the 64-bit filing system is 1, so the returned value is
    a minimum of 256 for a 64-bit file and less for a 32-bit file.
    \param nFid An integer file handle (1-8)
    \return (Major_Rev * 256)+Minor_Rev
    */
    MATINT_API int S64GetVersion(const int nFid);

    //! Get the offset to the next block (units of DBSize) the file would allocate
    /*!
    When a file is opened for reading, the physical file size should be the same or up to
    DBSize+DLSize less than this. If it is more, the file has extra information written on
    the end and needs fixing... it was probably interrupted during writing. This does not
    included buffered data that is not committed to disk during writing.
    \param nFid An integer file handle (1-8)
    \return The (estimated) file size in bytes. Should be a multiple of DBSize.
    */
    MATINT_API long long S64FileSize(const int nFid);

    //! Get the maximum time of any item in the file
    /*!
    The maximum time written to disk is saved in the file head. If this is not set the maximum
    time is found by scanning all the channels for the maximum time saved in any channel. We
    can also force all channels to be scanned with the bReadChans argument.
    \param nFid An integer file handle (1-8)
    \return The maximum time in the file in ticks. If the file is empty, this can be
    returned as -1.
    */
    MATINT_API long long S64MaxTime(const int nFid);


    //! 
    /*!

    */
    MATINT_API int S64TimeDate(const int nFid, long long *pTDGet, const long long *pTDSet, int iMode);

    //! 
    /*!

    */
    MATINT_API int S64AppID(const int nFid, int *pTDGet, const int *pTDSet, int iMode);

    MATINT_API int S64GetExtraData(const int nFid, void *pData, unsigned int nBytes, unsigned int nOffset);

    MATINT_API int S64SetExtraData(const int nFid, const void *pData, unsigned int nBytes, unsigned int nOffset);

    /*Get/Set Channel properties for the file*/

    //! Get the type of a channel
    /*!
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \return An integer corresponding the to type of the channel
    0 - No channel
    1 - Waveform
    2 - Event Fall
    3 - Event Rise
    4 - Event Both
    5 - Marker
    6 - Wavemark
    7 - Realmark
    8 - Textmark
    9 - Realwave
    */
    MATINT_API int S64ChanType(const int nFid, const int nChan);

    //! Get the channel divide rate
    /*!
    The channel divide rate is the number of of clocks ticks per sample interval
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \return The sample interval in file clock ticks. 
    */
    MATINT_API long long S64ChanDivide(const int nFid, const int nChan);

    //! Get the "ideal" rate for a channel
    /*!
    The ideal rate for a waveform channel should be set to be the desired rate, which will
    often be the reciprocal of the divide (scaled for seconds). The ChanDivide() routine
    will give you the actual rate. For an event-based channel, this value represents the
    expected mean sustained event rate - often used for scaling buffer allocations and also
    for guessing the y axis range for rate-related displays.

    This function gets the ideal waveform sample rate for Adc and
    RealWave channels and the expected average event rate for all other channel. 
    \param nFid An integer file handle (1-8)
    \param nChan The channel number
    \return     The previous value of the ideal rate or 0 if the channel does not exist.
    */
    MATINT_API double S64GetIdealRate(const int nFid, const int nChan);

    //! Set the "ideal" rate for a channel
    /*!
    This function sets the ideal waveform sample rate for Adc and
    RealWave channels and the expected average event rate for all other channel. 
    \param nFid An integer file handle (1-8)
    \param nChan The channel number
    \param dRate If this is >= 0.0 the value will be set. If negative, no change is made.
    \return     The previous value of the ideal rate or 0 if the channel does not exist.
    */
    MATINT_API double S64SetIdealRate(const int nFid, const int nChan, const double dRate);

    //! Get the comment text associated with a channel
    /*!
    Used to get a channel comment. This function copies data from one char buffer to another. It is
    *CRITICAL* that the buffers are the right size or *MATLAB WILL CRASH*. What you have to do call this
    function with nGetSize set to a negative number. The function wont copy any data, but will return
    the size of the required buffer. Then, in MATLAB, use this number to create a buffer of the correct
    size and call the function again with nGetSize set to a non-negative number.

    See S64SetChanComment() for more information.
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \param Comment A string holding the comment to set.
    \param nGetSize A flag to determine the behaviour of the function. If ngetSize < 0, the function
    returns the length of the comment, but does not get the comment. If nGetSize >= 0 the function
    gets the comment
    \return If ngetSize < 0, the length of the comment of a negative error message (the buffer needs 
    to be this number +1)
    If nGetSize >= 0, 0 or a negative code.
    */
    MATINT_API int S64GetChanComment(const int nFid, const int nChan, char *Comment, const int nGetSize);

    //! Set a channel comment
    /*!
    \param nFid An integer file handle (1-8)
    \param nChan The channel number
    \param Comment A string holding the comment to set.
    \return 0 or a negative code.
    */
    MATINT_API int S64SetChanComment(const int nFid, const int nChan, const char *Comment);

    //! Get the channel title
    /*!
    Used to get a channel title. This function copies data from one char buffer to another. It is
    *CRITICAL* that the buffers are the right size or *MATLAB WILL CRASH*. What you have to do call this
    function with nGetSize set to a negative number. The function wont copy any data, but will return
    the size of the required buffer. Then, in MATLAB, use this number to create a buffer of the correct
    size and call the function again with nGetSize set to a non-negative number.

    See S64SetChanTitle() for more information.
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \param Title A string holding the title to set.
    \param nGetSize A flag to determine the behaviour of the function. If ngetSize < 0, the function
    returns the length of the title, but does not get the comment. If nGetSize >= 0 the function
    gets the comment
    \return If ngetSize < 0, the length of the title of a negative error message (the buffer needs 
    to be this number +1)
    If nGetSize >= 0, 0 or a negative code.
    */
    MATINT_API int S64GetChanTitle(const int nFid, const int nChan, char *Title, const int nGetSize);

    //! Set a channel title
    /*!
    \param nFid An integer file handle (1-8)
    \param nChan The channel number
    \param Title A string holding the title to set.
    \return 0 or a negative code.
    */
    MATINT_API int S64SetChanTitle(const int nFid, const int nChan, const char *Title);

    //! Get the channel scale
    /*!
    See S64SetChanScale() for details.
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (can be any channel type as long as it exists)
    \param dScale The channel scale is returned here.
    \return S64_OK (0) or NO_CHANNEL if chan does not exist.
    */
    MATINT_API int S64GetChanScale(const int nFid, const int nChan, double *dScale);

    //! Set the channel scale
    /*!
    Channel scales are used to translate between integer representations of values and
    real units. They are used for Adc, RealWave and AdcMark channels and could find uses
    in other situations in the future. When a channel is expressed in user units:

    user units = integer * scale / 6553.6 + offset

    For a RealWave channel, where the channel is already in user units, the scale and
    offset values tell us how to convert the channel back into integers:

    integer = (user units - offset)*6553.6/scale

    The factor of 6553.6 comes about because if you have a 16-bit ADC that spans 10 V, a
    scale value of 1.0 converts between Volts and ADC values. Please do not worry about
    this value... just use the equations and all will be well.

    See also: S64GetChanScale(), S64SetChanOffset(), S64GetChanOffset()
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (can be any channel type as long as it exists)
    \param dScale The new scale value to apply. Any finite value except 0 is acceptable.
    However, using values close to the largest or smallest allowed fp values
    is stupid and risks over or underflows.
    \return S64_OK (0) or NO_CHANNEL if chan does not exist.
    */
    MATINT_API int S64SetChanScale(const int nFid, const int nChan, const double dScale);

    //! Get the channel offset
    /*!
    See S64SetChanScale() for details. See also S64SetChanoffset().
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (can be any channel type as long as it exists)
    \param dOffset The new channel offset in user units.
    \return S64_OK (0) or NO_CHANNEL if chan does not exist.
    */
    MATINT_API int S64GetChanOffset(const int nFid, const int nChan, double *dOffset);

    //! Set the channel offset
    /*!
    See S64SetChanScale() for details. See also S64GetChanScale().
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (can be any channel type as long as it exists)
    \param dOffset The new channel offset in user units.
    \return S64_OK (0) or NO_CHANNEL if chan does not exist.
    */
    MATINT_API int S64SetChanOffset(const int nFid, const int nChan, const double dOffset);

    //! Get the channel units
    /*!
    Used to get a channel units. This function copies data from one char buffer to another. It is
    *CRITICAL* that the buffers are the right size or *MATLAB WILL CRASH*. What you have to do call this
    function with nGetSize set to a negative number. The function wont copy any data, but will return
    the size of the required buffer. Then, in MATLAB, use this number to create a buffer of the correct
    size and call the function again with nGetSize set to a non-negative number.

    See S64SetChanUnits() for more information.
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \param Units A string holding the unit to set.
    \param nGetSize A flag to determine the behaviour of the function. If ngetSize < 0, the function
    returns the length of the unit, but does not get the comment. If nGetSize >= 0 the function
    gets the comment
    \return If ngetSize < 0, the length of the unit of a negative error message (the buffer needs 
    to be this number +1)
    If nGetSize >= 0, 0 or a negative code.
    */
    MATINT_API int S64GetChanUnits(const int nFid, const int nChan, char *Units, const int nGetSize);

    //! Set a channel title
    /*!
    \param nFid An integer file handle (1-8)
    \param nChan The channel number
    \param Units A string holding the units to set.
    \return 0 or a negative code.
    */    
    MATINT_API int S64SetChanUnits(const int nFid, const int nChan, const char *Units);

    //! Get the time of the last item on the channel
    /*!
    Returns the time in ticks of the last item held in the channel or -1 if no data exists.
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \returns The last time in the channel or -1 if no data (or file!) or NO_CHANNEL.
    */
    MATINT_API long long S64ChanMaxTime(const int nFid, const int nChan);

    //! Get the time of N items before a given time
    /*!
    This is easy to describe for an event channel. For a waveform channel, it is the
    latest time at which a read of N items would end before the nominated time. This
    may mean that a waveform read would not get N items. The problem with waveforms is
    due to gaps.
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \param tFrom The time before which we are to search.
    \param tTo The earliest time to search.
    \param n     The number of points backwards to search. If you need more than 32-bits worth or
    search you will have to loop this... and ask yourself why!
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return -1 if no such item found, a negative error code (NO_CHANNEL) or the time in ticks.
    */
    MATINT_API long long S64PrevNTime(const int nFid, const int nChan, long long tFrom, long long tTo,
        const int n, const int nMask, const int nAsWave);

    //! Delete a channel
    /*!
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    */
    MATINT_API int S64ChanDelete(const int nFid, const int nChan);

    //! (Attempt to) Undelete a channel
    /*!
    This command deals with undeleting a channel. To find if you can undelete a channel you ask for
    the old channel type, which will be ChanOff if you cannot delete it as in use or not deleted.
    Otherwise it will be the type of the channel that you can undelete.
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \return 0 if all OK, otherwise a negative error code.
    */
    MATINT_API int S64ChanUndelete(const int nFid, const int nChan);
    MATINT_API int S64GetChanYRange(const int nFid, const int nChan, double *dLow, double *dHigh);
    MATINT_API int S64SetChanYRange(const int nFid, const int nChan,
        const double dLow, const double dHigh);

    //! Get the byte size of the repeating object held by the channel
    /*!
    Each channel holds repeating objects, all the same size. Note that all extended
    marker types have sizes that are rounded up to a multiple of 8 bytes. This routine
    is a convenient way to get the object sizes (though you could calculate them). The
    non-Marker derived classes have fixed sizes, but these are also reported:

    channel type | size 
    -------------| ------
    Adc          |   2
    RealWave     |   4
    EventRise    |   8
    EventFall    |   8
    EventBoth    |  16
    Marker       |  16
    RealMark     |  16 + nRows * nCols * 4
    TextMark     |  16 + nRows
    AdcMark      |  16 + nRows * nCols * 2

    \param nFid An integer file handle (1-8)
    \param nChan The channel number.
    \return NO_CHANNEL or the size in bytes. The sizes of all extended marker types are rounded
    up to the next multiple of 8 bytes.
    */
    MATINT_API int S64ItemSize(const int nFid, const int nChan);

    /*Event Channels*/

    //! Create an EventFall or an EventRise channel
    /*!
    A basic event is a 64-bit (long long) time. EventBoth data is slightly different (it is essentially
    a marker) and should be created with the S64SetLevelChan() command.

    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1) for an unused channel
    \param dRate The expected maximum sustained (over several seconds) event rate.
    \param iType The type of event channel 1 = EventFall, 2 = EventRise, 3 = EventBoth
    We accept iType = 3 for the sake of completeness, but in this case the call will be
    re-routed to a S64SetLevelChan() call
    \return 0 or a negative code.
    */
    MATINT_API int S64SetEventChan(const int nFid, const int nChan, const double dRate, const int iType);

    //! Write event data to event channel
    /*!
    Write EventFall or EventRise data to an event channel. Use WriteLevels() for an EventBoth
    channel. The data must be after all data written previously to the channel. The data must
    be in ascending time order (or the file will be corrupt).

    \param nFid An integer file handle (1-8)
    \param nChan  An EventFall or EventRise channel.
    \param pData A buffer of data (int64) to be written
    \param nCount The number of event times to write
    \return S64_OK (0) or a negative error code (READ_ONLY, NO_CHANNEL, BAD_WRITE)
    */
    MATINT_API int S64WriteEvents(const int nFid, const int nChan, const long long *pData, const int nCount);

    //! Read event times from any suitable channel
    /*!
    You can read event times from any event channel, marker channel or extended marker channel.

    See also: SetEventChan(), WriteEvents(), ReadLevels()
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \param pData A buffer to accept the data of size at least nMax.
    \param nMax  The maximun number of events to read
    \param tFrom Earliest time we are interested in
    \param tTo One tick beyond the times we are interested in (up to but not including)
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return The number of events read or a negative error code
    */
    MATINT_API int S64ReadEvents(const int nFid, const int nChan,  long long *pData, int nMax,
        const long long tFrom, const long long tTo, const int nMask);

    /*Marker channels*/

    //! Create a Marker or a EventBoth channel
    /*!
    Markers are 64-bit (long long) time stamps plus marker codes. In the 64-bit son library
    EventBoth data is stored as a Marker, not as an event, using code 0 in the first dcoe to mean low
    and any other code to mean high.

    See also: SetLevelChan(), SetMarkerChan(), WriteMarkers(), ReadMarkers(), SetBuffering()
    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1) for an unused channel
    \param dRate    The expected maximum sustained (over several seconds) event rate. 
    \param nkind     The type of the channel, 5 = marker channel, 4 = level.
    \return 0 if no error is detected, otherwise a negative error code.
    */
    MATINT_API int S64SetMarkerChan(const int nFid, const int nChan, const double dRate, const int nkind);

    //! Write Marker data to a Marker or EventBoth channel.
    /*!
    You can also use UseS64 WriteLevels() for an EventBoth channel. The data must be after
    all data written previously to the channel. The data must be in ascending time order
    (or the file will be corrupt).
    \param nFid An integer file handle (1-8)
    \param nChan  An EventBoth or marker channel.
    \param pData A buffer of data to be written. Times must be in ascending order and the
    data must occur after any existing data in the channel.
    \param count The number of event times to write
    \return 0 or a negative error code
    */
    MATINT_API int S64WriteMarkers(const int nFid, const int nChan, const S64Marker* pData, const int count);

    //! Read Marker data from a Marker or extended Marker channel
    /*!
    You can read marker data from any channel that contains markers.
    \param nFid An integer file handle (1-8)
    \param nChan A Marker, EventBoth or extended marker channel.
    \param pData A buffer to receive the marker data.
    \param nMax  The maximum number of Markers to return.
    \param tFrom The first time to include in the search for markers.
    \param tUpto The first time not to include. Returned data will span the time range
    tRom up to but not including tUpto.
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return The number of markers read or a negative error code.
    */
    MATINT_API int S64ReadMarkers(const int nFid, const int nChan, S64Marker* pData, const int nMax,
        const long long tFrom, const long long tUpto, const int nMask);

    //! Modify a marker data item already in the channel
    /*!
    This is exclusively used to change marker codes of data already written to
    a Marker-based channel (but NOT the time of an item).
    \param nFid An integer file handle (1-8)
    \param nChan     An EventBoth or marker channel.
    \param t         The time of the existng marker to be modified
    \param pM        The marker object to replace the current object. We stromgly recommend that
                     the time is the same as t (and may decde to ignore the new time).
    */
    MATINT_API int S64EditMarker(const int nFid, const int nChan, long long t, const S64Marker* pM);

    /*Marker Masks*/

    //! Read a marker mask to a 256-by-4 matrix
    /*!
    \param nMask     The integer mask handle.
    \param iCode     The mode of the mask (0 = AND, 1 = OR)
    \param nMode     An array of ints to contain the 256-by-4 matrix data
    */
    MATINT_API int S64GetMaskCodes(const int nMask, int* iCode, int* nMode);

    //! Set the codes in a marker mask
    /*!
    \param nMask     The integer mask handle.
    \param nCode     An integer 256-by-4 matrix. If the i-jth entry is non-zero,
                     then item i in layer j of the mask is included in the mask
    */
    MATINT_API int S64SetMaskCodes(const int nMask, const int* nCode);

    //! Set the mode of a marker mask
    /*!
    \param nMask     The integer mask handle.
    \param nMode     The mode of the mask (0 = AND, 1 = OR)
    */
    MATINT_API int S64SetMaskMode(const int nMask, const int nMode);

    //! Get the mode of a marker mask
    /*!
    \param nMask     The integer mask handle.
    \param nMode    The mode of the mask (0 = AND, 1 = OR)
    */
    MATINT_API int S64GetMaskMode(const int nMask);

    //! Set the column select of a marker mask
    /*!
    \param nMask     The integer mask handle.
    \param nCol      The new column select value, or -1 to reset
    */
    MATINT_API int S64SetMaskCol(const int nMask, const int nCol);

    //! Get the column select of a marker mask
    /*!
    \param nMask     The integer mask handle.
    */
    MATINT_API int S64GetMaskCol(const int nMask);

    //! Resets marker mask by removing all codes from all layers
    /*!
    \param nMask     The integer mask handle.
    */
    MATINT_API int S64ResetMask(const int nMask);

    //! Resets all marker mask by removing all codes from all layers in all filters
    /*!
    */
    MATINT_API int S64ResetAllMasks();

    /*Level Channels*/

    //! Create an EventBoth channel
    /*!
    In the 64-bit SON library we store EventBoth data as Marker data using code 0 of the
    first marker code for low and any other code for high. This command sets the initial
    level to low, equivalent to S64SetInitLevel(nFid, nChan, 0).
    \param nFid      An integer file handle (1-8)
    \param nChan     An unused channel in the file.
    \param dRate     The expected maximum sustained event rate (used to allocate buffers)
    \return 0 or a negative error code.
    */
    MATINT_API int S64SetLevelChan(const int nFid, const int nChan, const double dRate);

    //! Set the initial level of an EventBoth channel
    /*!
    This command only has any effect if used after creating a channel and before any data
    is written with S64WriteLevels(). It sets the level of the input before the first change
    and so defines the direction of the first change.

    \param nFid      An integer file handle (1-8)
    \param nChan     An existing EventBoth channel
    \param nLevel    Set =0 if the initial level is low, !=0 if it is high. This sets the
                     input level at time -1 ticks and this level is assumed to continue until
                     data is written to the channel.
    \return 0 or a negative error code.
    */
    MATINT_API int S64SetInitLevel(const int nFid, const int nChan, const int nLevel);

    //! Write EventBoth data assuming that each time toggles the level
    /*!
    The written data must lie after all preceding data and the first
    time is assumed to have the opposite level to the last written level. If no data has been
    written, the first level is assumed to toggle the level set by SetInitLevel(). 

    \param nFid      An integer file handle (1-8)
    \param nChan     An existing EventBoth channel
    \param pData     A buffer of event times that are in ascending time order and that occur
                     after any data already written to the file.
    \param nCount    The number of events to add to the file.
    \return      0 or a negative error code (READ_ONLY, NO_CHANNEL, BAD_WRITE)
    */
    MATINT_API int S64WriteLevels(const int nFid, const int nChan, const long long *pData, const int nCount);

    //! Read alternating levels from an EventBoth channel
    /*!
    This is provided for compatibility with the 32-bit SON library. It assumes that the channel
    holds alternating high and low levels and reads it back as events, not Markers.

    \param nFid      An integer file handle (1-8)
    \param nChan     An existing EventBoth channel
    \param pData     A buffer to hold any returned times.
    \param nMax      The maximum number of times to return.
    \param tFrom     The first time to include in the search for levels.
    \param tTo       The first time not to include. Returned data will span the time range
                     tRom up to but not including tUpto.
    \param nLevel    The level of the first returned point (1 if high, 0 if low) or if no
                     points are in the time range, the level at time tFrom.
    \return The numbe rof points read or a negative error code.
    */
    MATINT_API int S64ReadLevels(const int nFid, const int nChan, long long *pData, int nMax,
        const long long tFrom, const long long tTo, int* nLevel);

    /*Extended Markers*/

    //! Create a TextMark channel
    /*!
    A TextMark is a Marker plus a zero terminated 8-bit character string. Data is written and read
    using S64Write1TextMark() and S64Read1TextMark(). Note that as the size of this object is not fixed,
    you should use S64GetExtMarkInfo() to find the size. This is very important when writing or reading
    Extended markers as passing in the wrong sized buffers will crash MATLAB

    You can create this channel type with the S64SetExtMarkChan() command, if you wish.

    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1) for an unused channel
    \param dRate The expected sustained maximum event rate.
    \param nMax  The number of 8-bit bytes (characters) attached to each TextMark, including the zero byte
    to terminate the string. In the 64-bit library setting this to a multiple of 8 makes use
    of all the allocated space.
    \return 0 or a negative error code.
    */
    MATINT_API int S64SetTextMarkChan(const int nFid, const int nChan, double dRate, const int nMax);

    //! Create a TextMark, RealMark or AdcMark channel
    /*!
    This command creates and extended marker channel, that is a channel of Markers with attached
    other data. This data is arranged as a grid of rows and columns of data. The data is 16-bit
    integers for AdcMark, 32-bit float values for RealMark and 8-bit characters for TextMark. 

    The 64-bit SON library rounds up the size of an extended marker to a multiple of 8 bytes. The
    32-bit SON library rounded up the size to a multiple of 4 bytes. Data is written and read
    using S64Write1XXXMark() and S64Read1XXXMark(). Note that as the size of this object is not fixed,
    You can use ItemSize() to find the size to move points on by and GetExtMarkInfo() where
    the rows is equivalent to the nMax parameter.

    S64SetTextMarkChan() is a special case of the command to create a TextMark channel with a
    single column.

    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1) for an unused channel
    \param dRate The expected sustained maximum event rate.
    \param nType  1 = AdcMark (short items), 2 = RealMark (float items) or 3 = TextMark (chars)
    \param nRows The number of rows (at least 1) of attached items.
    \param nCols The number of columns (at least 1) of attached items.
    \param tDiv  If > 0 this means that the attached data (taken by row) has a time axis.
    This is usually used for AdcMark data, but it could also be used for RealMark
    channels. The data in each row is assumed to start at the Marker time and be
    spaced by tDiv. ChanDivide() returns this value. tDiv is in units of the file time base.
    \return 0 or a negative error code.
    */
    MATINT_API int S64SetExtMarkChan(const int nFid, const int nChan, double dRate, const int nType, 
        const int nRows, const int nCols, const long long tDiv);

    //! Get extended marker data information
    /*!
    Use this command to read the number of rows and columns.

    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1)
    \param nRows Either nullptr or points at a variable to hold the number of rows
    \param nCols Either nullptr or points at a variable to hold the number of columns
    \return The number of pre-alignment points or a negative error code.
    */
    MATINT_API int S64GetExtMarkInfo(const int nFid, const int nChan, int* nRows, int* nCols);

    //! Writes a single Textmarker to a Textmark channel
    /*!

    \param nFid An integer file handle (1-8)
    \param nChan    The channel number (starts at 1)
    \param pData    The basic marker data for the textmarker
    \param text     the string
    \param nSize    the length of the the string
    \return 0 or a negative error code.
    */
    MATINT_API int S64Write1TextMark(const int nFid, const int nChan, const S64Marker* pData,
        const char* text, const int nSize);

    //! Writes a single Realmarker to a realmark channel
    /*!

    \param nFid An integer file handle (1-8)
    \param nChan    The channel number (starts at 1)
    \param pData    The basic marker data for the realmarker
    \param real     the realmarker data
    \param nSize    the the number of reals
    \return 0 or a negative error code.
    */
    MATINT_API int S64Write1RealMark(const int nFid, const int nChan, const S64Marker* pData,
        const float* real, const int nSize);

    //! Writes a single Wavemarker to a wavemark channel
    /*!

    \param nFid An integer file handle (1-8)
    \param nChan    The channel number (starts at 1)
    \param pData    The basic marker data for the wavemarker
    \param wave     the wavemarker data
    \param nSize    the the number of shorts
    \return 0 or a negative error code.
    */
    MATINT_API int S64Write1WaveMark(const int nFid, const int nChan, const S64Marker* pData,
        const short* wave, const int nSize);


    //! Reads a single Textmarker from a Textmark channel
    /*!
    Reads the first textmarker in the interval tFrom to tUp (both in clock ticks)

    \param nFid An integer file handle (1-8)
    \param nChan    The channel number (starts at 1)
    \param pData    The basic marker data for the textmarker
    \param text     the string
    \param tFrom    the start time
    \param tUpto    the end time
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return 0 or a negative error code.
    */
    MATINT_API int S64Read1TextMark(const int nFid, const int nChan, S64Marker* pData, char* text,
        const long long tFrom, const long long tUpto, const int nMask);

    //! Reads a single Realmarker from a realmark channel
    /*!
    Reads the first Realmarker in the interval tFrom to tUp (both in clock ticks)

    \param nFid An integer file handle (1-8)
    \param nChan    The channel number (starts at 1)
    \param pData    The basic marker data for the realmarker
    \param real     the realmarker data
    \param tFrom    the start time
    \param tUpto    the end time
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return 0 or a negative error code.
    */
    MATINT_API int S64Read1RealMark(const int nFid, const int nChan, S64Marker* pData, float* real,
        const long long tFrom, const long long tUpto, const int nMask);

    //! Reads a single Wavemarker from a wavemark channel
    /*!
    Reads the first Wavemarker in the interval tFrom to tUp (both in clock ticks)

    \param nFid An integer file handle (1-8)
    \param nChan    The channel number (starts at 1)
    \param pData    The basic marker data for the wavemarker
    \param real     the wavemarker data
    \param tFrom    the start time
    \param tUpto    the end time
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return 0 or a negative error code.
    */
    MATINT_API int S64Read1WaveMark(const int nFid, const int nChan, S64Marker* pData, short* real,
        const long long tFrom, const long long tUpto, const int nMask);

    /*Wave Channels*/

    //! Create a waveform or realwave channel using short or float data items
    /*!
    The use of short data (16-bit signed integers) is convenient and compact (waveform data usually acounts for
    the bulk of the disk space used for data files. It also matches a common ADC specification. However, user
    data is more conveniently used in user units, so waveform channels have an asscoiated scale and offset value
    that is used to convert between integer units and real units.

    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1) for an unused channel
    \param tDiv The spacing of the channel in timebase units (S64SetTimeBase(), S64GetTimeBase())
    \param nType Either 1 for Adc or 9 for RealWave, the channel type.
    \param dRate The desired sample rate in Hz. It can happen that due to the choice of the file time base, the
    lDvd value is only an approximation to the desired rate. This stores the desired rate for
    information purposes, see IdealRate().
    \return 0 or a negative error code.
    */
    MATINT_API int S64SetWaveChan(const int nFid, const int nChan, const long long tDiv, const int nType,
        const double dRate);

    //! Write waveform data as shorts to an Adc channel
    /*!
    Unlike event-based channel types where all data must be written after any data already present in a channel,
    we allow you to overwrite previously-written Adc data.

    In normal use, all waveform data for a channel will be aligned at time that are the first data time plus an
    integer value times the lDvd for the channel. We do not prevent you writing data where the alignment is not
    the same after a gap, but programs like Spike2 will have subtle problems if you do this.


    \param nFid An integer file handle (1-8)
    \param nChan The channel number (starts at 1) This must be an Adc (not a waveform) channel.
    \param pData The buffer of data to write.
    \param count The number of values to write.
    \param tFrom This is the time of the first item in the buffer to write. If this is after the last written data,
    new values are appended to the file. If it is at or before the last time written, old data is
    replaced.
    \return The next time to write to or a negative error code.
    */
    MATINT_API long long S64WriteWaveS(const int nFid, const int nChan, const short* pData,
        const int count, const long long tFrom);

    //! Write waveform data as floats to a RealWave channel
    /*!
    Unlike event-based channel types where all data must be written after any data already present in a channel,
    we allow you to overwrite previously-written RealWave data.

    In normal use, all waveform data for a channel will be aligned at time that are the first data time plus an
    integer value times the lDvd for the channel. We do not prevent you writing data where the alignment is not
    the same after a gap, but programs like Spike2 will have subtle problems if you do this.

    \param nFid An integer file handle (1-8)
    \param nChan The waveform channel to write to. This must be a RealWave (not an ADC) channel.
    \param pData The buffer of data to write.
    \param count The number of values to write.
    \param tFrom This is the time of the first item in the buffer to write. If this is after the last written data,
    new values are appended to the file. If it is at or before the last time written, old data is
    replaced.
    \return The next time to write to or a negative error code.
    */
    MATINT_API long long S64WriteWaveF(const int nFid, const int nChan,  const float* pData,
        const int count, const long long tFrom);

    //! Write waveform data as floats to a RealWave channel
    /*!
    Unlike event-based channel types where all data must be written after any data already present in a channel,
    we allow you to overwrite previously-written RealWave data.

    In normal use, all waveform data for a channel will be aligned at time that are the first data time plus an
    integer value times the lDvd for the channel. We do not prevent you writing data where the alignment is not
    the same after a gap, but programs like Spike2 will have subtle problems if you do this.

    Like S64ReadWave64() this will not give you any greater accuracy since we immeadiately cast
    to floats. This function should only be used if for some reason pData absolutely has to be buffer of
    doubles.

    \param nFid An integer file handle (1-8)
    \param nChan The waveform channel to write to. This must be a RealWave (not an ADC) channel.
    \param pData The buffer of data to write.
    \param count The number of values to write.
    \param tFrom This is the time of the first item in the buffer to write. If this is after the last written data,
    new values are appended to the file. If it is at or before the last time written, old data is
    replaced.
    \return The next time to write to or a negative error code.
    */
    MATINT_API long long S64WriteWave64(const int nFid, const int nChan, const double* pData,
        const int count, const long long tFrom);

    //! Read an Adc, RealWave or an AdcMark channel as shorts
    /*!
    If you read a RealWave channel, the values are converted to integers, the values are converted to 
    shorts using the scale and offset set for the channel.
    The read operation stops at a gap in the data. A gap is defined as an interval between data points 
    that is not the lDvd value defined for the channel.

    \param nFid An integer file handle (1-8)
    \param nChan The waveform channel to read. This can be either Adc or RealWave. 
    \param pData The buffer to receive the read data.
    \param nMax The maximum number of values to read.
    \param tFrom This and tUpto define a time range in which to locate the data to read. The first data point will
    be at or after tFrom.
    \param tUpto Not data returned will be at or after this time.
    \param tFirst If any data is returned this is set to the time of the first item.
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return The number of values returned or a negative error code.
    */
    MATINT_API int S64ReadWaveS(const int nFid, const int nChan, short* pData, const int nMax,
        const long long tFrom, const long long tUpto, long long* tFirst, const int nMask);

    //! Read an Adc, RealWave or an AdcMark channel as floats
    /*!
    If you read an Adc or AdcMark channel, the values are converted to floats using the scale and offset
    set for the channel.
    The read operation stops at a gap in the data. A gap is defined as an interval between data points 
    that is not the lDvd value defined for the channel.

    \param nFid An integer file handle (1-8)
    \param nChan The waveform channel to read. This can be either Adc or RealWave. 
    \param pData The buffer to receive the read data.
    \param nMax The maximum number of values to read.
    \param tFrom This and tUpto define a time range in which to locate the data to read. The first data point will
    be at or after tFrom.
    \param tUpto Not data returned will be at or after this time.
    \param tFirst If any data is returned this is set to the time of the first item.
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return The number of values returned or a negative error code.
    */
    MATINT_API int S64ReadWaveF(const int nFid, const int nChan, float* pData, const int nMax,
        const long long tFrom, const long long tUpto, long long* tFirst, const int nMask);

    //! Read an Adc, RealWave or an AdcMark channel as doubles
    /*!
    The values are converted from shorts or floats to doubles using the scale and offset set 
    for the channel. Note this will not give you any extra accuracy or any other benefit. The only
    effect this will have is using more memory, this function should only be used if for some reason
    you absolutely need the resulting data set to be doubles.

    The read operation stops at a gap in the data. A gap is defined as an interval between 
    data points that is not the lDvd value defined for the channel.

    \param nFid An integer file handle (1-8)
    \param nChan The waveform channel to read. This can be either Adc or RealWave. If the type is RealWave we convert
    values to short using the channel scale and offset. This may result in truncated values.
    \param pData The buffer to receive the read data.
    \param nMax The maximum number of values to read.
    \param tFrom This and tUpto define a time range in which to locate the data to read. The first data point will
    be at or after tFrom.
    \param tUpto Not data returned will be at or after this time.
    \param tFirst If any data is returned this is set to the time of the first item.
    \param nMask Integer filter handle Used with AdcMark channels to filter the data.
    \return The number of values returned or a negative error code.
    */
    MATINT_API int S64ReadWave64(const int nFid, const int nChan, double* pData, const int nMax,
        const long long tFrom, const long long tUpto, long long* tFirst, const int nMask);

#ifdef __cplusplus
}
#endif

#endif
