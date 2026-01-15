"""
Valuations Library Core Bindings
"""

from __future__ import annotations

import collections.abc
import typing

import numpy
import numpy.typing

__all__: list[str] = [
    "Actual360",
    "Actual365Fixed",
    "Actual365NoLeap",
    "ActualActual",
    "AfSabr",
    "American",
    "AntonovSabr",
    "Argentina",
    "Australia",
    "Bermudan",
    "BespokeCalendar",
    "Brazil",
    "Business252",
    "CMS",
    "CTree",
    "CalendarType",
    "Call",
    "Canada",
    "China",
    "CmsSpread",
    "CzechRepublic",
    "DayCounterType",
    "Denmark",
    "European",
    "ExerciseStyle",
    "Finland",
    "Germany",
    "HHW",
    "HWPDE",
    "Heston",
    "HongKong",
    "HullWhite",
    "Hungary",
    "Iceland",
    "India",
    "Indonesia",
    "Italy",
    "Japan",
    "LocalVol",
    "LogNormalBasket",
    "Mexico",
    "NewZealand",
    "NonRecombining",
    "Norway",
    "OneDayCounter",
    "OptionType",
    "Poland",
    "Put",
    "QuantoedCMS",
    "QuantoedCmsSpread",
    "Recombining",
    "Russia",
    "Sabr",
    "SabrPDE",
    "SaudiArabia",
    "SchobelZhu",
    "ShortRate1FPDE",
    "ShortRate2FPDE",
    "SimpleDayCounter",
    "Singapore",
    "SkewMC",
    "Slovakia",
    "SouthAfrica",
    "SouthKorea",
    "Swaption",
    "Sweden",
    "Switzerland",
    "TARGET",
    "Taiwan",
    "Termstructure",
    "Thirty360",
    "TreeType",
    "Turkey",
    "Ukraine",
    "UnitedKingdom",
    "black_formula_call",
    "black_formula_call_vega",
    "black_scholes_call",
    "black_scholes_call_vega",
    "cdf_normal",
    "implied_vol",
    "pdf_normal",
]

class AfSabr(SabrPDE):
    def __init__(
        self,
        alpha: typing.SupportsFloat,
        beta: typing.SupportsFloat,
        nu: typing.SupportsFloat,
        rho: typing.SupportsFloat,
        shift: typing.SupportsFloat,
        maturity: typing.SupportsFloat,
        F: typing.SupportsFloat,
        sizeX: typing.SupportsInt,
        sizeT: typing.SupportsInt,
        nd: typing.SupportsFloat,
    ) -> None: ...
    def getShift(self) -> float: ...

class AntonovSabr(SabrPDE):
    def __init__(
        self,
        alpha: typing.SupportsFloat,
        beta: typing.SupportsFloat,
        nu: typing.SupportsFloat,
        rho: typing.SupportsFloat,
        maturity: typing.SupportsFloat,
        F: typing.SupportsFloat,
        sizeX: typing.SupportsInt,
        sizeT: typing.SupportsInt,
        nd: typing.SupportsFloat,
    ) -> None: ...

class CMS:
    def __init__(
        self,
        expirySR: typing.SupportsFloat,
        tenorSR: typing.SupportsFloat,
        forwardSR: typing.SupportsFloat,
        annuitySR: typing.SupportsFloat,
        payCMS: typing.SupportsFloat,
        discountCMS: typing.SupportsFloat,
        beta: typing.SupportsFloat = 0.85,
        alpha: typing.SupportsFloat = 0.5,
        nu: typing.SupportsFloat = 0.25,
        rho: typing.SupportsFloat = -0.75,
    ) -> None: ...
    def fairValue(
        self, strike: typing.SupportsFloat, callORput: str = "call"
    ) -> float: ...
    def getATMvol(self) -> float: ...
    def getDiscountCMS(self) -> float: ...
    def getForward(self) -> float: ...
    def getImpliedVol(self, strike: typing.SupportsFloat) -> float: ...
    def getMaturity(self) -> float: ...
    @property
    def alpha(self) -> float: ...
    @property
    def nu(self) -> float: ...
    @property
    def rho(self) -> float: ...

class CTree:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        S: typing.SupportsFloat,
        T: collections.abc.Sequence[typing.SupportsFloat],
        F: collections.abc.Sequence[typing.SupportsFloat],
        IV: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...
    def calculateBinomial(
        self,
        strike: typing.SupportsFloat,
        Maturity: typing.SupportsFloat,
        Nnodes: typing.SupportsInt,
        style: ExerciseStyle,
        pay: OptionType,
        tree: TreeType,
    ) -> float: ...
    def calculateTrinomial(
        self,
        strike: typing.SupportsFloat,
        Maturity: typing.SupportsFloat,
        Nnodes: typing.SupportsInt,
        style: ExerciseStyle,
        pay: OptionType,
        tree: TreeType,
    ) -> float: ...

class CalendarType:
    """
    Members:

      UnitedKingdom

      Argentina

      Australia

      BespokeCalendar

      Brazil

      Canada

      China

      CzechRepublic

      Denmark

      Finland

      Germany

      HongKong

      Hungary

      Iceland

      India

      Indonesia

      Italy

      Japan

      Ukraine

      Mexico

      NewZealand

      Norway

      Poland

      Russia

      SaudiArabia

      Singapore

      Slovakia

      SouthAfrica

      SouthKorea

      Sweden

      Switzerland

      Taiwan

      TARGET

      Turkey
    """

    Argentina: typing.ClassVar[CalendarType]  # value = <CalendarType.Argentina: 2>
    Australia: typing.ClassVar[CalendarType]  # value = <CalendarType.Australia: 3>
    BespokeCalendar: typing.ClassVar[
        CalendarType
    ]  # value = <CalendarType.BespokeCalendar: 4>
    Brazil: typing.ClassVar[CalendarType]  # value = <CalendarType.Brazil: 5>
    Canada: typing.ClassVar[CalendarType]  # value = <CalendarType.Canada: 6>
    China: typing.ClassVar[CalendarType]  # value = <CalendarType.China: 7>
    CzechRepublic: typing.ClassVar[
        CalendarType
    ]  # value = <CalendarType.CzechRepublic: 8>
    Denmark: typing.ClassVar[CalendarType]  # value = <CalendarType.Denmark: 9>
    Finland: typing.ClassVar[CalendarType]  # value = <CalendarType.Finland: 10>
    Germany: typing.ClassVar[CalendarType]  # value = <CalendarType.Germany: 11>
    HongKong: typing.ClassVar[CalendarType]  # value = <CalendarType.HongKong: 12>
    Hungary: typing.ClassVar[CalendarType]  # value = <CalendarType.Hungary: 13>
    Iceland: typing.ClassVar[CalendarType]  # value = <CalendarType.Iceland: 14>
    India: typing.ClassVar[CalendarType]  # value = <CalendarType.India: 15>
    Indonesia: typing.ClassVar[CalendarType]  # value = <CalendarType.Indonesia: 16>
    Italy: typing.ClassVar[CalendarType]  # value = <CalendarType.Italy: 17>
    Japan: typing.ClassVar[CalendarType]  # value = <CalendarType.Japan: 18>
    Mexico: typing.ClassVar[CalendarType]  # value = <CalendarType.Mexico: 20>
    NewZealand: typing.ClassVar[CalendarType]  # value = <CalendarType.NewZealand: 21>
    Norway: typing.ClassVar[CalendarType]  # value = <CalendarType.Norway: 22>
    Poland: typing.ClassVar[CalendarType]  # value = <CalendarType.Poland: 23>
    Russia: typing.ClassVar[CalendarType]  # value = <CalendarType.Russia: 24>
    SaudiArabia: typing.ClassVar[CalendarType]  # value = <CalendarType.SaudiArabia: 25>
    Singapore: typing.ClassVar[CalendarType]  # value = <CalendarType.Singapore: 26>
    Slovakia: typing.ClassVar[CalendarType]  # value = <CalendarType.Slovakia: 27>
    SouthAfrica: typing.ClassVar[CalendarType]  # value = <CalendarType.SouthAfrica: 28>
    SouthKorea: typing.ClassVar[CalendarType]  # value = <CalendarType.SouthKorea: 29>
    Sweden: typing.ClassVar[CalendarType]  # value = <CalendarType.Sweden: 30>
    Switzerland: typing.ClassVar[CalendarType]  # value = <CalendarType.Switzerland: 31>
    TARGET: typing.ClassVar[CalendarType]  # value = <CalendarType.TARGET: 33>
    Taiwan: typing.ClassVar[CalendarType]  # value = <CalendarType.Taiwan: 32>
    Turkey: typing.ClassVar[CalendarType]  # value = <CalendarType.Turkey: 34>
    Ukraine: typing.ClassVar[CalendarType]  # value = <CalendarType.Ukraine: 19>
    UnitedKingdom: typing.ClassVar[
        CalendarType
    ]  # value = <CalendarType.UnitedKingdom: 1>
    __members__: typing.ClassVar[
        dict[str, CalendarType]
    ]  # value = {'UnitedKingdom': <CalendarType.UnitedKingdom: 1>, 'Argentina': <CalendarType.Argentina: 2>, 'Australia': <CalendarType.Australia: 3>, 'BespokeCalendar': <CalendarType.BespokeCalendar: 4>, 'Brazil': <CalendarType.Brazil: 5>, 'Canada': <CalendarType.Canada: 6>, 'China': <CalendarType.China: 7>, 'CzechRepublic': <CalendarType.CzechRepublic: 8>, 'Denmark': <CalendarType.Denmark: 9>, 'Finland': <CalendarType.Finland: 10>, 'Germany': <CalendarType.Germany: 11>, 'HongKong': <CalendarType.HongKong: 12>, 'Hungary': <CalendarType.Hungary: 13>, 'Iceland': <CalendarType.Iceland: 14>, 'India': <CalendarType.India: 15>, 'Indonesia': <CalendarType.Indonesia: 16>, 'Italy': <CalendarType.Italy: 17>, 'Japan': <CalendarType.Japan: 18>, 'Ukraine': <CalendarType.Ukraine: 19>, 'Mexico': <CalendarType.Mexico: 20>, 'NewZealand': <CalendarType.NewZealand: 21>, 'Norway': <CalendarType.Norway: 22>, 'Poland': <CalendarType.Poland: 23>, 'Russia': <CalendarType.Russia: 24>, 'SaudiArabia': <CalendarType.SaudiArabia: 25>, 'Singapore': <CalendarType.Singapore: 26>, 'Slovakia': <CalendarType.Slovakia: 27>, 'SouthAfrica': <CalendarType.SouthAfrica: 28>, 'SouthKorea': <CalendarType.SouthKorea: 29>, 'Sweden': <CalendarType.Sweden: 30>, 'Switzerland': <CalendarType.Switzerland: 31>, 'Taiwan': <CalendarType.Taiwan: 32>, 'TARGET': <CalendarType.TARGET: 33>, 'Turkey': <CalendarType.Turkey: 34>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class CmsSpread:
    def __init__(
        self,
        expiry1: typing.SupportsFloat,
        tenor1: typing.SupportsFloat,
        fwd1: typing.SupportsFloat,
        annuity1: typing.SupportsFloat,
        pay1: typing.SupportsFloat,
        disc1: typing.SupportsFloat,
        beta1: typing.SupportsFloat,
        strikes1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        quotes1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        type1: str,
        expiry2: typing.SupportsFloat,
        tenor2: typing.SupportsFloat,
        fwd2: typing.SupportsFloat,
        annuity2: typing.SupportsFloat,
        pay2: typing.SupportsFloat,
        disc2: typing.SupportsFloat,
        beta2: typing.SupportsFloat,
        strikes2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        quotes2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        type2: str,
        corr: typing.SupportsFloat,
    ) -> None: ...
    def simulationCMSs(self) -> list[float]: ...
    def spreadOption(
        self, K: typing.SupportsFloat, a: typing.SupportsFloat, b: typing.SupportsFloat
    ) -> float: ...

class DayCounterType:
    """
    Members:

      Actual360

      Actual365Fixed

      ActualActual

      Actual365NoLeap

      Business252

      OneDayCounter

      SimpleDayCounter

      Thirty360
    """

    Actual360: typing.ClassVar[DayCounterType]  # value = <DayCounterType.Actual360: 1>
    Actual365Fixed: typing.ClassVar[
        DayCounterType
    ]  # value = <DayCounterType.Actual365Fixed: 2>
    Actual365NoLeap: typing.ClassVar[
        DayCounterType
    ]  # value = <DayCounterType.Actual365NoLeap: 4>
    ActualActual: typing.ClassVar[
        DayCounterType
    ]  # value = <DayCounterType.ActualActual: 3>
    Business252: typing.ClassVar[
        DayCounterType
    ]  # value = <DayCounterType.Business252: 5>
    OneDayCounter: typing.ClassVar[
        DayCounterType
    ]  # value = <DayCounterType.OneDayCounter: 6>
    SimpleDayCounter: typing.ClassVar[
        DayCounterType
    ]  # value = <DayCounterType.SimpleDayCounter: 7>
    Thirty360: typing.ClassVar[DayCounterType]  # value = <DayCounterType.Thirty360: 8>
    __members__: typing.ClassVar[
        dict[str, DayCounterType]
    ]  # value = {'Actual360': <DayCounterType.Actual360: 1>, 'Actual365Fixed': <DayCounterType.Actual365Fixed: 2>, 'ActualActual': <DayCounterType.ActualActual: 3>, 'Actual365NoLeap': <DayCounterType.Actual365NoLeap: 4>, 'Business252': <DayCounterType.Business252: 5>, 'OneDayCounter': <DayCounterType.OneDayCounter: 6>, 'SimpleDayCounter': <DayCounterType.SimpleDayCounter: 7>, 'Thirty360': <DayCounterType.Thirty360: 8>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class ExerciseStyle:
    """
    Members:

      American

      European

      Bermudan
    """

    American: typing.ClassVar[ExerciseStyle]  # value = <ExerciseStyle.American: 0>
    Bermudan: typing.ClassVar[ExerciseStyle]  # value = <ExerciseStyle.Bermudan: 2>
    European: typing.ClassVar[ExerciseStyle]  # value = <ExerciseStyle.European: 1>
    __members__: typing.ClassVar[
        dict[str, ExerciseStyle]
    ]  # value = {'American': <ExerciseStyle.American: 0>, 'European': <ExerciseStyle.European: 1>, 'Bermudan': <ExerciseStyle.Bermudan: 2>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HHW:
    def HHWPrice(
        self, maturity: typing.SupportsFloat, strike: typing.SupportsFloat
    ) -> float: ...
    def __init__(
        self,
        s0: typing.SupportsFloat,
        v0: typing.SupportsFloat,
        r0: typing.SupportsFloat,
        kappa: typing.SupportsFloat,
        eta: typing.SupportsFloat,
        rho: typing.SupportsFloat,
        sigma1: typing.SupportsFloat,
        sigma2: typing.SupportsFloat,
        a: typing.SupportsFloat,
    ) -> None: ...

class HWPDE:
    @typing.overload
    def __init__(
        self,
        R0: typing.SupportsFloat,
        kappa: typing.SupportsFloat,
        timeSigmas: collections.abc.Sequence[typing.SupportsFloat],
        sigmas: collections.abc.Sequence[typing.SupportsFloat],
        timeThetas: collections.abc.Sequence[typing.SupportsFloat],
        thetas: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        kappa: typing.SupportsFloat,
        timeSigmas: collections.abc.Sequence[typing.SupportsFloat],
        sigmas: collections.abc.Sequence[typing.SupportsFloat],
        timeDF: collections.abc.Sequence[typing.SupportsFloat],
        DF: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...
    def pricingBermudan(
        self,
        Expiry: typing.SupportsFloat,
        Tenor: typing.SupportsFloat,
        Exercises: collections.abc.Sequence[typing.SupportsFloat],
        Strike: typing.SupportsFloat,
        PayFrequency: typing.SupportsFloat,
    ) -> float: ...
    def pricingCallableSwap(
        self,
        Expiry: typing.SupportsFloat,
        Tenor: typing.SupportsFloat,
        Exercises: collections.abc.Sequence[typing.SupportsFloat],
        Coupon: typing.SupportsFloat,
        Strike: typing.SupportsFloat,
        PayFrequency: typing.SupportsFloat,
        type: str,
    ) -> float: ...
    def pricingSwaption(
        self,
        Expiry: typing.SupportsFloat,
        Tenor: typing.SupportsFloat,
        Strike: typing.SupportsFloat,
        PayFrequency: typing.SupportsFloat,
    ) -> float: ...

class Heston:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        spot: typing.SupportsFloat,
        var0: typing.SupportsFloat,
        kappa: typing.SupportsFloat,
        theta: typing.SupportsFloat,
        xi: typing.SupportsFloat,
        rho: typing.SupportsFloat,
        seed: typing.SupportsInt = 42,
    ) -> None: ...
    def hestonPrice(
        self,
        maturity: typing.SupportsFloat,
        forward: typing.SupportsFloat,
        strike: typing.SupportsFloat,
        optType: str = "call",
    ) -> float: ...
    def simulationHeston(
        self,
        times: collections.abc.Sequence[typing.SupportsFloat],
        forwards: collections.abc.Sequence[typing.SupportsFloat],
    ) -> list[float]: ...
    @property
    def kappa(self) -> float: ...
    @property
    def rho(self) -> float: ...
    @property
    def theta(self) -> float: ...
    @property
    def var0(self) -> float: ...
    @property
    def xi(self) -> float: ...

class HullWhite:
    def __init__(
        self,
        kappa: typing.SupportsFloat,
        timeSigmas: collections.abc.Sequence[typing.SupportsFloat],
        sigmas: collections.abc.Sequence[typing.SupportsFloat],
        timeDFs: collections.abc.Sequence[typing.SupportsFloat],
        DFs: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...
    def getParameterKappa(self) -> float: ...
    def getSigmas(self) -> list[float]: ...
    def getTimeSigmas(self) -> list[float]: ...
    def optionBond(
        self,
        Expiry: typing.SupportsFloat,
        Maturity: typing.SupportsFloat,
        Strike: typing.SupportsFloat,
        callORput: str,
    ) -> float: ...
    def simulationHW(
        self, times: collections.abc.Sequence[typing.SupportsFloat]
    ) -> list[float]: ...
    def swaption(
        self,
        Expiry: typing.SupportsFloat,
        Tenor: typing.SupportsFloat,
        Strike: typing.SupportsFloat,
        PayFrequency: typing.SupportsFloat = 0.5,
    ) -> float: ...

class LocalVol:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, sabrModels: collections.abc.Sequence[Sabr]) -> None: ...
    def callPDE(
        self,
        maturity: typing.SupportsFloat,
        strike: typing.SupportsFloat,
        N: typing.SupportsInt = 100,
    ) -> float: ...
    def density(
        self, maturity: typing.SupportsFloat, Nt: typing.SupportsInt
    ) -> list[float]: ...
    def putPDE(
        self,
        maturity: typing.SupportsFloat,
        strike: typing.SupportsFloat,
        N: typing.SupportsInt = 100,
    ) -> float: ...
    @property
    def spot(self) -> float: ...
    @spot.setter
    def spot(self, arg1: typing.SupportsFloat) -> None: ...

class LogNormalBasket:
    def __init__(
        self,
        Spot: collections.abc.Sequence[typing.SupportsFloat],
        Strike: collections.abc.Sequence[typing.SupportsFloat],
        Maturities: collections.abc.Sequence[typing.SupportsFloat],
        Forwards: collections.abc.Sequence[
            collections.abc.Sequence[typing.SupportsFloat]
        ],
        IV: collections.abc.Sequence[collections.abc.Sequence[typing.SupportsFloat]],
        correlation: collections.abc.Sequence[
            collections.abc.Sequence[typing.SupportsFloat]
        ],
    ) -> None: ...
    def get_nassets(self) -> int: ...
    def simulate_basket(
        self, schedule: collections.abc.Sequence[typing.SupportsFloat]
    ) -> list[list[float]]: ...
    def simulate_basketWR(
        self, schedule: collections.abc.Sequence[typing.SupportsFloat]
    ) -> list[float]: ...

class OptionType:
    """
    Members:

      Call

      Put
    """

    Call: typing.ClassVar[OptionType]  # value = <OptionType.Call: 0>
    Put: typing.ClassVar[OptionType]  # value = <OptionType.Put: 1>
    __members__: typing.ClassVar[
        dict[str, OptionType]
    ]  # value = {'Call': <OptionType.Call: 0>, 'Put': <OptionType.Put: 1>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class QuantoedCMS:
    def __init__(
        self,
        expiry: typing.SupportsFloat,
        tenor: typing.SupportsFloat,
        fwd: typing.SupportsFloat,
        annuity: typing.SupportsFloat,
        pay: typing.SupportsFloat,
        disc: typing.SupportsFloat,
        corFX: typing.SupportsFloat,
        atmVolFX: typing.SupportsFloat,
        beta: typing.SupportsFloat,
        strikes: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        quotes: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        type: str = "premium",
    ) -> None: ...
    def fairValue(
        self, strike: typing.SupportsFloat, callORput: str = "call"
    ) -> float: ...
    def getForward(self) -> float: ...
    def simulation(self, corrRN: typing.SupportsFloat) -> float: ...

class QuantoedCmsSpread:
    def __init__(
        self,
        expiry1: typing.SupportsFloat,
        tenor1: typing.SupportsFloat,
        fwd1: typing.SupportsFloat,
        annuity1: typing.SupportsFloat,
        pay1: typing.SupportsFloat,
        disc1: typing.SupportsFloat,
        corFX1: typing.SupportsFloat,
        atmVolFX1: typing.SupportsFloat,
        beta1: typing.SupportsFloat,
        strikes1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        quotes1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        type1: str,
        expiry2: typing.SupportsFloat,
        tenor2: typing.SupportsFloat,
        fwd2: typing.SupportsFloat,
        annuity2: typing.SupportsFloat,
        pay2: typing.SupportsFloat,
        disc2: typing.SupportsFloat,
        corFX2: typing.SupportsFloat,
        atmVolFX2: typing.SupportsFloat,
        beta2: typing.SupportsFloat,
        strikes2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        quotes2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"],
        type2: str,
        corr: typing.SupportsFloat,
    ) -> None: ...
    def simulationQuantoedCMSs(
        self, cr1: typing.SupportsFloat, cr2: typing.SupportsFloat
    ) -> list[float]: ...

class Sabr:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        maturity: typing.SupportsFloat,
        forward: typing.SupportsFloat,
        beta: typing.SupportsFloat = 0.85,
        alpha: typing.SupportsFloat = 0.5,
        nu: typing.SupportsFloat = 0.25,
        rho: typing.SupportsFloat = -0.75,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        maturity: typing.SupportsFloat,
        forward: typing.SupportsFloat,
        beta: typing.SupportsFloat,
        alpha: typing.SupportsFloat,
        nu: typing.SupportsFloat,
        rho: typing.SupportsFloat,
        shift: typing.SupportsFloat,
    ) -> None: ...
    def impliedVol(self, strike: typing.SupportsFloat) -> float: ...
    def localVol(self, spot: typing.SupportsFloat) -> float: ...
    def normalVol(self, K: typing.SupportsFloat) -> float: ...
    def premiumBachelier(
        self, strike: typing.SupportsFloat, callORput: str = "call"
    ) -> float: ...
    def premiumBlackScholes(
        self, strike: typing.SupportsFloat, callORput: str = "call"
    ) -> float: ...
    @property
    def alpha(self) -> float: ...
    @alpha.setter
    def alpha(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def beta(self) -> float: ...
    @beta.setter
    def beta(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def forward(self) -> float: ...
    @forward.setter
    def forward(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def maturity(self) -> float: ...
    @maturity.setter
    def maturity(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def nu(self) -> float: ...
    @nu.setter
    def nu(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def rho(self) -> float: ...
    @rho.setter
    def rho(self, arg1: typing.SupportsFloat) -> None: ...

class SabrPDE:
    def getAlpha(self) -> float: ...
    def getBeta(self) -> float: ...
    def getDensity(self) -> list[float]: ...
    def getFgrid(self) -> list[float]: ...
    def getNu(self) -> float: ...
    def getRho(self) -> float: ...

class SchobelZhu:
    def SchobelPrice(
        self,
        maturity: typing.SupportsFloat,
        forward: typing.SupportsFloat,
        strike: typing.SupportsFloat,
    ) -> float: ...
    def __init__(
        self,
        spot: typing.SupportsFloat,
        var0: typing.SupportsFloat,
        kappa: typing.SupportsFloat,
        theta: typing.SupportsFloat,
        xi: typing.SupportsFloat,
        rho: typing.SupportsFloat,
    ) -> None: ...
    def calibrator(
        self,
        maturitys: collections.abc.Sequence[typing.SupportsFloat],
        forwards: collections.abc.Sequence[typing.SupportsFloat],
        strikes: collections.abc.Sequence[typing.SupportsFloat],
        marketQuotes: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...
    def simulationSchobelZhu(
        self,
        times: collections.abc.Sequence[typing.SupportsFloat],
        forwards: collections.abc.Sequence[typing.SupportsFloat],
    ) -> list[float]: ...
    @property
    def kappa(self) -> float: ...
    @kappa.setter
    def kappa(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def rho(self) -> float: ...
    @rho.setter
    def rho(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def theta(self) -> float: ...
    @theta.setter
    def theta(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def var0(self) -> float: ...
    @var0.setter
    def var0(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def xi(self) -> float: ...
    @xi.setter
    def xi(self, arg1: typing.SupportsFloat) -> None: ...

class ShortRate1FPDE:
    def __init__(
        self,
        R0: typing.SupportsFloat,
        kappa: typing.SupportsFloat,
        alpha: typing.SupportsFloat,
        beta: typing.SupportsFloat,
        gamma: typing.SupportsFloat,
        timeSigmas: collections.abc.Sequence[typing.SupportsFloat],
        sigmas: collections.abc.Sequence[typing.SupportsFloat],
        timeThetas: collections.abc.Sequence[typing.SupportsFloat],
        thetas: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...
    def pricingSwaption(
        self,
        Expiry: typing.SupportsFloat,
        Tenor: typing.SupportsFloat,
        Strike: typing.SupportsFloat,
        PayFrequency: typing.SupportsFloat = 0.5,
    ) -> float: ...

class ShortRate2FPDE:
    def __init__(
        self,
        kappa1: typing.SupportsFloat,
        kappa2: typing.SupportsFloat,
        lambda_: typing.SupportsFloat,
        timeSigma1s: collections.abc.Sequence[typing.SupportsFloat],
        sigma1s: collections.abc.Sequence[typing.SupportsFloat],
        timeSigma2s: collections.abc.Sequence[typing.SupportsFloat],
        sigma2s: collections.abc.Sequence[typing.SupportsFloat],
        timeAlphas: collections.abc.Sequence[typing.SupportsFloat],
        alphas: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...
    def pricingSwaption(
        self,
        Expiry: typing.SupportsFloat,
        Tenor: typing.SupportsFloat,
        Strike: typing.SupportsFloat,
        PayFrequency: typing.SupportsFloat = 0.5,
    ) -> float: ...

class SkewMC:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, sabrModels: collections.abc.Sequence[Sabr]) -> None: ...
    def simulation(
        self,
        times: collections.abc.Sequence[typing.SupportsFloat],
        spot: typing.SupportsFloat,
        kappa: typing.SupportsFloat,
    ) -> list[float]: ...

class Swaption:
    def __init__(
        self,
        expiry: typing.SupportsFloat,
        tenor: typing.SupportsFloat,
        forward: typing.SupportsFloat,
        annuity: typing.SupportsFloat,
        beta: typing.SupportsFloat = 0.85,
        alpha: typing.SupportsFloat = 0.5,
        nu: typing.SupportsFloat = 0.25,
        rho: typing.SupportsFloat = -0.75,
    ) -> None: ...
    def getImpliedVol(self, strike: typing.SupportsFloat) -> float: ...
    def swapFairValue(self, strike: typing.SupportsFloat) -> float: ...
    def swaptionFairValue(
        self, strike: typing.SupportsFloat, callORput: str = "call"
    ) -> float: ...
    @property
    def alpha(self) -> float: ...
    @property
    def nu(self) -> float: ...
    @property
    def rho(self) -> float: ...

class Termstructure:
    @typing.overload
    def __init__(
        self,
        days: collections.abc.Sequence[typing.SupportsFloat],
        rate: collections.abc.Sequence[typing.SupportsFloat],
        calendar: CalendarType,
        daycount: DayCounterType,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        days: collections.abc.Sequence[typing.SupportsFloat],
        rate: collections.abc.Sequence[typing.SupportsFloat],
        qdays: collections.abc.Sequence[typing.SupportsFloat],
        divident: collections.abc.Sequence[typing.SupportsFloat],
        calendar: CalendarType,
        daycount: DayCounterType,
    ) -> None: ...
    def discount(self, date: typing.SupportsInt) -> float: ...
    def divident(
        self, date: typing.SupportsInt, tenor: typing.SupportsInt
    ) -> float: ...
    def rate(self, date: typing.SupportsInt, tenor: typing.SupportsInt) -> float: ...

class TreeType:
    """
    Members:

      Recombining

      NonRecombining
    """

    NonRecombining: typing.ClassVar[TreeType]  # value = <TreeType.NonRecombining: 1>
    Recombining: typing.ClassVar[TreeType]  # value = <TreeType.Recombining: 0>
    __members__: typing.ClassVar[
        dict[str, TreeType]
    ]  # value = {'Recombining': <TreeType.Recombining: 0>, 'NonRecombining': <TreeType.NonRecombining: 1>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

def black_formula_call(
    Forward: typing.SupportsFloat,
    Strike: typing.SupportsFloat,
    Vol: typing.SupportsFloat,
    Expiry: typing.SupportsFloat,
) -> float: ...
def black_formula_call_vega(
    Forward: typing.SupportsFloat,
    Strike: typing.SupportsFloat,
    Vol: typing.SupportsFloat,
    Expiry: typing.SupportsFloat,
) -> float: ...
def black_scholes_call(
    Spot: typing.SupportsFloat,
    Strike: typing.SupportsFloat,
    r: typing.SupportsFloat,
    d: typing.SupportsFloat,
    Vol: typing.SupportsFloat,
    Expiry: typing.SupportsFloat,
) -> float: ...
def black_scholes_call_vega(
    Spot: typing.SupportsFloat,
    Strike: typing.SupportsFloat,
    r: typing.SupportsFloat,
    d: typing.SupportsFloat,
    Vol: typing.SupportsFloat,
    Expiry: typing.SupportsFloat,
) -> float: ...
def cdf_normal(arg0: typing.SupportsFloat) -> float:
    """
    Cumulative distribution function for normal distribution
    """

def implied_vol(
    Maturity: typing.SupportsFloat,
    Forward: typing.SupportsFloat,
    Strike: typing.SupportsFloat,
    Price: typing.SupportsFloat,
) -> float:
    """
    Calculate implied volatility for a Black-Scholes price
    """

def pdf_normal(arg0: typing.SupportsFloat) -> float:
    """
    Probability density function for normal distribution
    """

Actual360: DayCounterType  # value = <DayCounterType.Actual360: 1>
Actual365Fixed: DayCounterType  # value = <DayCounterType.Actual365Fixed: 2>
Actual365NoLeap: DayCounterType  # value = <DayCounterType.Actual365NoLeap: 4>
ActualActual: DayCounterType  # value = <DayCounterType.ActualActual: 3>
American: ExerciseStyle  # value = <ExerciseStyle.American: 0>
Argentina: CalendarType  # value = <CalendarType.Argentina: 2>
Australia: CalendarType  # value = <CalendarType.Australia: 3>
Bermudan: ExerciseStyle  # value = <ExerciseStyle.Bermudan: 2>
BespokeCalendar: CalendarType  # value = <CalendarType.BespokeCalendar: 4>
Brazil: CalendarType  # value = <CalendarType.Brazil: 5>
Business252: DayCounterType  # value = <DayCounterType.Business252: 5>
Call: OptionType  # value = <OptionType.Call: 0>
Canada: CalendarType  # value = <CalendarType.Canada: 6>
China: CalendarType  # value = <CalendarType.China: 7>
CzechRepublic: CalendarType  # value = <CalendarType.CzechRepublic: 8>
Denmark: CalendarType  # value = <CalendarType.Denmark: 9>
European: ExerciseStyle  # value = <ExerciseStyle.European: 1>
Finland: CalendarType  # value = <CalendarType.Finland: 10>
Germany: CalendarType  # value = <CalendarType.Germany: 11>
HongKong: CalendarType  # value = <CalendarType.HongKong: 12>
Hungary: CalendarType  # value = <CalendarType.Hungary: 13>
Iceland: CalendarType  # value = <CalendarType.Iceland: 14>
India: CalendarType  # value = <CalendarType.India: 15>
Indonesia: CalendarType  # value = <CalendarType.Indonesia: 16>
Italy: CalendarType  # value = <CalendarType.Italy: 17>
Japan: CalendarType  # value = <CalendarType.Japan: 18>
Mexico: CalendarType  # value = <CalendarType.Mexico: 20>
NewZealand: CalendarType  # value = <CalendarType.NewZealand: 21>
NonRecombining: TreeType  # value = <TreeType.NonRecombining: 1>
Norway: CalendarType  # value = <CalendarType.Norway: 22>
OneDayCounter: DayCounterType  # value = <DayCounterType.OneDayCounter: 6>
Poland: CalendarType  # value = <CalendarType.Poland: 23>
Put: OptionType  # value = <OptionType.Put: 1>
Recombining: TreeType  # value = <TreeType.Recombining: 0>
Russia: CalendarType  # value = <CalendarType.Russia: 24>
SaudiArabia: CalendarType  # value = <CalendarType.SaudiArabia: 25>
SimpleDayCounter: DayCounterType  # value = <DayCounterType.SimpleDayCounter: 7>
Singapore: CalendarType  # value = <CalendarType.Singapore: 26>
Slovakia: CalendarType  # value = <CalendarType.Slovakia: 27>
SouthAfrica: CalendarType  # value = <CalendarType.SouthAfrica: 28>
SouthKorea: CalendarType  # value = <CalendarType.SouthKorea: 29>
Sweden: CalendarType  # value = <CalendarType.Sweden: 30>
Switzerland: CalendarType  # value = <CalendarType.Switzerland: 31>
TARGET: CalendarType  # value = <CalendarType.TARGET: 33>
Taiwan: CalendarType  # value = <CalendarType.Taiwan: 32>
Thirty360: DayCounterType  # value = <DayCounterType.Thirty360: 8>
Turkey: CalendarType  # value = <CalendarType.Turkey: 34>
Ukraine: CalendarType  # value = <CalendarType.Ukraine: 19>
UnitedKingdom: CalendarType  # value = <CalendarType.UnitedKingdom: 1>
