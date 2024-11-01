/*!
 * Paper.js v0.9.20 - The Swiss Army Knife of Vector Graphics Scripting.
 * http://paperjs.org/
 *
 * Copyright (c) 2011 - 2014, Juerg Lehni & Jonathan Puckey
 * http://scratchdisk.com/ & http://jonathanpuckey.com/
 *
 * Distributed under the MIT license. See LICENSE file for details.
 *
 * All rights reserved.
 *
 * Date: Mon Aug 25 14:21:13 2014 +0200
 *
 ***
 *
 * Straps.js - Class inheritance library with support for bean-style accessors
 *
 * Copyright (c) 2006 - 2013 Juerg Lehni
 * http://scratchdisk.com/
 *
 * Distributed under the MIT license.
 *
 ***
 *
 * Acorn.js
 * http://marijnhaverbeke.nl/acorn/
 *
 * Acorn is a tiny, fast JavaScript parser written in JavaScript,
 * created by Marijn Haverbeke and released under an MIT license.
 *
 */
var paper = new (function (t) {
    var e = new (function () {
        function n(t, n, i, r, a) {
            function o(s, o) {
                (o = o || ((o = u(n, s)) && (o.get ? o : o.value))),
                    "string" == typeof o &&
                        "#" === o[0] &&
                        (o = t[o.substring(1)] || o);
                var l,
                    d = "function" == typeof o,
                    f = o,
                    _ = a || d ? (o && o.get ? s in t : t[s]) : null;
                (a && _) ||
                    (d && _ && (o.base = _),
                    d &&
                        r !== !1 &&
                        (l = s.match(/^([gs]et|is)(([A-Z])(.*))$/)) &&
                        (h[l[3].toLowerCase() + l[4]] = l[2]),
                    (f &&
                        !d &&
                        f.get &&
                        "function" == typeof f.get &&
                        e.isPlainObject(f)) ||
                        (f = {
                            value: f,
                            writable: !0
                        }),
                    (
                        u(t, s) || {
                            configurable: !0
                        }
                    ).configurable &&
                        ((f.configurable = !0), (f.enumerable = i)),
                    c(t, s, f));
            }
            var h = {};
            if (n) {
                for (var l in n) n.hasOwnProperty(l) && !s.test(l) && o(l);
                for (var l in h) {
                    var d = h[l],
                        f = t["set" + d],
                        _ = t["get" + d] || (f && t["is" + d]);
                    !_ ||
                        (r !== !0 && 0 !== _.length) ||
                        o(l, {
                            get: _,
                            set: f
                        });
                }
            }
            return t;
        }
        function i(t, e, n) {
            return (
                t &&
                    ("length" in t &&
                    !t.getLength &&
                    "number" == typeof t.length
                        ? a
                        : o
                    ).call(t, e, (n = n || t)),
                n
            );
        }
        function r(t, e) {
            for (var n in e) e.hasOwnProperty(n) && (t[n] = e[n]);
            return t;
        }
        var s = /^(statics|enumerable|beans|preserve)$/,
            a =
                [].forEach ||
                function (t, e) {
                    for (var n = 0, i = this.length; i > n; n++)
                        t.call(e, this[n], n, this);
                },
            o = function (t, e) {
                for (var n in this)
                    this.hasOwnProperty(n) && t.call(e, this[n], n, this);
            },
            h =
                Object.create ||
                function (t) {
                    return {
                        __proto__: t
                    };
                },
            u =
                Object.getOwnPropertyDescriptor ||
                function (t, e) {
                    var n = t.__lookupGetter__ && t.__lookupGetter__(e);
                    return n
                        ? {
                              get: n,
                              set: t.__lookupSetter__(e),
                              enumerable: !0,
                              configurable: !0
                          }
                        : t.hasOwnProperty(e)
                        ? {
                              value: t[e],
                              enumerable: !0,
                              configurable: !0,
                              writable: !0
                          }
                        : null;
                },
            l =
                Object.defineProperty ||
                function (t, e, n) {
                    return (
                        (n.get || n.set) && t.__defineGetter__
                            ? (n.get && t.__defineGetter__(e, n.get),
                              n.set && t.__defineSetter__(e, n.set))
                            : (t[e] = n.value),
                        t
                    );
                },
            c = function (t, e, n) {
                return delete t[e], l(t, e, n);
            };
        return n(
            function () {
                for (var t = 0, e = arguments.length; e > t; t++)
                    r(this, arguments[t]);
            },
            {
                inject: function (t) {
                    if (t) {
                        var e = t.statics === !0 ? t : t.statics,
                            i = t.beans,
                            r = t.preserve;
                        e !== t && n(this.prototype, t, t.enumerable, i, r),
                            n(this, e, !0, i, r);
                    }
                    for (var s = 1, a = arguments.length; a > s; s++)
                        this.inject(arguments[s]);
                    return this;
                },
                extend: function () {
                    for (
                        var t, e = this, i = 0, r = arguments.length;
                        r > i && !(t = arguments[i].initialize);
                        i++
                    );
                    return (
                        (t =
                            t ||
                            function () {
                                e.apply(this, arguments);
                            }),
                        (t.prototype = h(this.prototype)),
                        (t.base = e),
                        c(t.prototype, "constructor", {
                            value: t,
                            writable: !0,
                            configurable: !0
                        }),
                        n(t, this, !0),
                        arguments.length ? this.inject.apply(t, arguments) : t
                    );
                }
            },
            !0
        ).inject({
            inject: function () {
                for (var t = 0, e = arguments.length; e > t; t++) {
                    var i = arguments[t];
                    i && n(this, i, i.enumerable, i.beans, i.preserve);
                }
                return this;
            },
            extend: function () {
                var t = h(this);
                return t.inject.apply(t, arguments);
            },
            each: function (t, e) {
                return i(this, t, e);
            },
            set: function (t) {
                return r(this, t);
            },
            clone: function () {
                return new this.constructor(this);
            },
            statics: {
                each: i,
                create: h,
                define: c,
                describe: u,
                set: r,
                clone: function (t) {
                    return r(new t.constructor(), t);
                },
                isPlainObject: function (t) {
                    var n = null != t && t.constructor;
                    return (
                        n && (n === Object || n === e || "Object" === n.name)
                    );
                },
                pick: function () {
                    for (var e = 0, n = arguments.length; n > e; e++)
                        if (arguments[e] !== t) return arguments[e];
                }
            }
        });
    })();
    "undefined" != typeof module && (module.exports = e),
        Array.isArray ||
            (Array.isArray = function (t) {
                return "[object Array]" === Object.prototype.toString.call(t);
            }),
        document.head ||
            (document.head = document.getElementsByTagName("head")[0]),
        e.inject({
            toString: function () {
                return null != this._id
                    ? (this._class || "Object") +
                          (this._name
                              ? " '" + this._name + "'"
                              : " @" + this._id)
                    : "{ " +
                          e
                              .each(
                                  this,
                                  function (t, e) {
                                      if (!/^_/.test(e)) {
                                          var n = typeof t;
                                          this.push(
                                              e +
                                                  ": " +
                                                  ("number" === n
                                                      ? a.instance.number(t)
                                                      : "string" === n
                                                      ? "'" + t + "'"
                                                      : t)
                                          );
                                      }
                                  },
                                  []
                              )
                              .join(", ") +
                          " }";
            },
            exportJSON: function (t) {
                return e.exportJSON(this, t);
            },
            toJSON: function () {
                return e.serialize(this);
            },
            _set: function (n, i, r) {
                if (n && (r || e.isPlainObject(n))) {
                    var s = n._filtering || n;
                    for (var a in s)
                        if (a in this && s.hasOwnProperty(a) && (!i || !i[a])) {
                            var o = n[a];
                            o !== t && (this[a] = o);
                        }
                    return !0;
                }
            },
            statics: {
                exports: {
                    enumerable: !0
                },
                extend: function ie() {
                    var t = ie.base.apply(this, arguments),
                        n = t.prototype._class;
                    return n && !e.exports[n] && (e.exports[n] = t), t;
                },
                equals: function (t, n) {
                    function i(t, e) {
                        for (var n in t)
                            if (t.hasOwnProperty(n) && !e.hasOwnProperty(n))
                                return !1;
                        return !0;
                    }
                    if (t === n) return !0;
                    if (t && t.equals) return t.equals(n);
                    if (n && n.equals) return n.equals(t);
                    if (Array.isArray(t) && Array.isArray(n)) {
                        if (t.length !== n.length) return !1;
                        for (var r = 0, s = t.length; s > r; r++)
                            if (!e.equals(t[r], n[r])) return !1;
                        return !0;
                    }
                    if (
                        t &&
                        "object" == typeof t &&
                        n &&
                        "object" == typeof n
                    ) {
                        if (!i(t, n) || !i(n, t)) return !1;
                        for (var r in t)
                            if (t.hasOwnProperty(r) && !e.equals(t[r], n[r]))
                                return !1;
                        return !0;
                    }
                    return !1;
                },
                read: function (n, i, r, s) {
                    if (this === e) {
                        var a = this.peek(n, i);
                        return n.__index++, a;
                    }
                    var o = this.prototype,
                        h = o._readIndex,
                        u = i || (h && n.__index) || 0;
                    s || (s = n.length - u);
                    var l = n[u];
                    return l instanceof this ||
                        (r && r.readNull && null == l && 1 >= s)
                        ? (h && (n.__index = u + 1),
                          l && r && r.clone ? l.clone() : l)
                        : ((l = e.create(this.prototype)),
                          h && (l.__read = !0),
                          (l =
                              l.initialize.apply(
                                  l,
                                  u > 0 || s < n.length
                                      ? Array.prototype.slice.call(n, u, u + s)
                                      : n
                              ) || l),
                          h && ((n.__index = u + l.__read), (l.__read = t)),
                          l);
                },
                peek: function (t, e) {
                    return t[(t.__index = e || t.__index || 0)];
                },
                remain: function (t) {
                    return t.length - (t.__index || 0);
                },
                readAll: function (t, e, n) {
                    for (var i, r = [], s = e || 0, a = t.length; a > s; s++)
                        r.push(
                            Array.isArray((i = t[s]))
                                ? this.read(i, 0, n)
                                : this.read(t, s, n, 1)
                        );
                    return r;
                },
                readNamed: function (n, i, r, s, a) {
                    var o = this.getNamed(n, i),
                        h = o !== t;
                    if (h) {
                        var u = n._filtered;
                        u ||
                            ((u = n._filtered = e.create(n[0])),
                            (u._filtering = n[0])),
                            (u[i] = t);
                    }
                    return this.read(h ? [o] : n, r, s, a);
                },
                getNamed: function (n, i) {
                    var r = n[0];
                    return (
                        n._hasObject === t &&
                            (n._hasObject =
                                1 === n.length && e.isPlainObject(r)),
                        n._hasObject ? (i ? r[i] : n._filtered || r) : t
                    );
                },
                hasNamed: function (t, e) {
                    return !!this.getNamed(t, e);
                },
                isPlainValue: function (t, e) {
                    return (
                        this.isPlainObject(t) ||
                        Array.isArray(t) ||
                        (e && "string" == typeof t)
                    );
                },
                serialize: function (t, n, i, r) {
                    n = n || {};
                    var s,
                        o = !r;
                    if (
                        (o &&
                            ((n.formatter = new a(n.precision)),
                            (r = {
                                length: 0,
                                definitions: {},
                                references: {},
                                add: function (t, e) {
                                    var n = "#" + t._id,
                                        i = this.references[n];
                                    if (!i) {
                                        this.length++;
                                        var r = e.call(t),
                                            s = t._class;
                                        s && r[0] !== s && r.unshift(s),
                                            (this.definitions[n] = r),
                                            (i = this.references[n] = [n]);
                                    }
                                    return i;
                                }
                            })),
                        t && t._serialize)
                    ) {
                        s = t._serialize(n, r);
                        var h = t._class;
                        !h || i || s._compact || s[0] === h || s.unshift(h);
                    } else if (Array.isArray(t)) {
                        s = [];
                        for (var u = 0, l = t.length; l > u; u++)
                            s[u] = e.serialize(t[u], n, i, r);
                        i && (s._compact = !0);
                    } else if (e.isPlainObject(t)) {
                        s = {};
                        for (var u in t)
                            t.hasOwnProperty(u) &&
                                (s[u] = e.serialize(t[u], n, i, r));
                    } else
                        s =
                            "number" == typeof t
                                ? n.formatter.number(t, n.precision)
                                : t;
                    return o && r.length > 0
                        ? [["dictionary", r.definitions], s]
                        : s;
                },
                deserialize: function (t, n, i) {
                    var r = t;
                    if (((i = i || {}), Array.isArray(t))) {
                        var s = t[0],
                            a = "dictionary" === s;
                        if (!a) {
                            if (i.dictionary && 1 == t.length && /^#/.test(s))
                                return i.dictionary[s];
                            s = e.exports[s];
                        }
                        r = [];
                        for (var o = s ? 1 : 0, h = t.length; h > o; o++)
                            r.push(e.deserialize(t[o], n, i));
                        if (a) i.dictionary = r[0];
                        else if (s) {
                            var u = r;
                            n
                                ? (r = n(s, u))
                                : ((r = e.create(s.prototype)), s.apply(r, u));
                        }
                    } else if (e.isPlainObject(t)) {
                        r = {};
                        for (var l in t) r[l] = e.deserialize(t[l], n, i);
                    }
                    return r;
                },
                exportJSON: function (t, n) {
                    var i = e.serialize(t, n);
                    return n && n.asString === !1 ? i : JSON.stringify(i);
                },
                importJSON: function (t, n) {
                    return e.deserialize(
                        "string" == typeof t ? JSON.parse(t) : t,
                        function (t, i) {
                            var r =
                                    n && n.constructor === t
                                        ? n
                                        : e.create(t.prototype),
                                s = r === n;
                            if (
                                1 === i.length &&
                                r instanceof y &&
                                (s || !(r instanceof x))
                            ) {
                                var a = i[0];
                                e.isPlainObject(a) && (a.insert = !1);
                            }
                            return t.apply(r, i), s && (n = null), r;
                        }
                    );
                },
                splice: function (e, n, i, r) {
                    var s = n && n.length,
                        a = i === t;
                    (i = a ? e.length : i), i > e.length && (i = e.length);
                    for (var o = 0; s > o; o++) n[o]._index = i + o;
                    if (a) return e.push.apply(e, n), [];
                    var h = [i, r];
                    n && h.push.apply(h, n);
                    for (
                        var u = e.splice.apply(e, h), o = 0, l = u.length;
                        l > o;
                        o++
                    )
                        u[o]._index = t;
                    for (var o = i + s, l = e.length; l > o; o++)
                        e[o]._index = o;
                    return u;
                },
                capitalize: function (t) {
                    return t.replace(/\b[a-z]/g, function (t) {
                        return t.toUpperCase();
                    });
                },
                camelize: function (t) {
                    return t.replace(/-(.)/g, function (t, e) {
                        return e.toUpperCase();
                    });
                },
                hyphenate: function (t) {
                    return t.replace(/([a-z])([A-Z])/g, "$1-$2").toLowerCase();
                }
            }
        });
    var n = {
            attach: function (n, i) {
                if ("string" != typeof n)
                    return (
                        e.each(
                            n,
                            function (t, e) {
                                this.attach(e, t);
                            },
                            this
                        ),
                        t
                    );
                var r = this._eventTypes[n];
                if (r) {
                    var s = (this._handlers = this._handlers || {});
                    (s = s[n] = s[n] || []),
                        -1 == s.indexOf(i) &&
                            (s.push(i),
                            r.install &&
                                1 == s.length &&
                                r.install.call(this, n));
                }
            },
            detach: function (n, i) {
                if ("string" != typeof n)
                    return (
                        e.each(
                            n,
                            function (t, e) {
                                this.detach(e, t);
                            },
                            this
                        ),
                        t
                    );
                var r,
                    s = this._eventTypes[n],
                    a = this._handlers && this._handlers[n];
                s &&
                    a &&
                    (!i || (-1 != (r = a.indexOf(i)) && 1 == a.length)
                        ? (s.uninstall && s.uninstall.call(this, n),
                          delete this._handlers[n])
                        : -1 != r && a.splice(r, 1));
            },
            once: function (t, e) {
                this.attach(t, function () {
                    e.apply(this, arguments), this.detach(t, e);
                });
            },
            fire: function (t, e) {
                var n = this._handlers && this._handlers[t];
                if (!n) return !1;
                for (
                    var i = [].slice.call(arguments, 1),
                        r = this,
                        s = 0,
                        a = n.length;
                    a > s;
                    s++
                )
                    if (n[s].apply(r, i) === !1 && e && e.stop) {
                        e.stop();
                        break;
                    }
                return !0;
            },
            responds: function (t) {
                return !(!this._handlers || !this._handlers[t]);
            },
            on: "#attach",
            off: "#detach",
            trigger: "#fire",
            _installEvents: function (t) {
                var e = this._handlers,
                    n = t ? "install" : "uninstall";
                for (var i in e)
                    if (e[i].length > 0) {
                        var r = this._eventTypes[i],
                            s = r[n];
                        s && s.call(this, i);
                    }
            },
            statics: {
                inject: function re() {
                    for (var t = 0, n = arguments.length; n > t; t++) {
                        var i = arguments[t],
                            r = i._events;
                        if (r) {
                            var s = {};
                            e.each(r, function (t, n) {
                                var r = "string" == typeof t,
                                    a = r ? t : n,
                                    o = e.capitalize(a),
                                    h = a.substring(2).toLowerCase();
                                (s[h] = r ? {} : t),
                                    (a = "_" + a),
                                    (i["get" + o] = function () {
                                        return this[a];
                                    }),
                                    (i["set" + o] = function (t) {
                                        var e = this[a];
                                        e && this.detach(h, e),
                                            t && this.attach(h, t),
                                            (this[a] = t);
                                    });
                            }),
                                (i._eventTypes = s);
                        }
                        re.base.call(this, i);
                    }
                    return this;
                }
            }
        },
        r = e.extend({
            _class: "PaperScope",
            initialize: function se() {
                if (
                    ((paper = this),
                    (this.settings = new e({
                        applyMatrix: !0,
                        handleSize: 4,
                        hitTolerance: 0
                    })),
                    (this.project = null),
                    (this.projects = []),
                    (this.tools = []),
                    (this.palettes = []),
                    (this._id = se._id++),
                    (se._scopes[this._id] = this),
                    !this.support)
                ) {
                    var t = Q.getContext(1, 1);
                    (se.prototype.support = {
                        nativeDash: "setLineDash" in t || "mozDash" in t,
                        nativeBlendModes: te.nativeModes
                    }),
                        Q.release(t);
                }
            },
            version: "0.9.20",
            getView: function () {
                return this.project && this.project.getView();
            },
            getPaper: function () {
                return this;
            },
            execute: function (t, e, n) {
                paper.PaperScript.execute(t, this, e, n), Z.updateFocus();
            },
            install: function (t) {
                var n = this;
                e.each(["project", "view", "tool"], function (i) {
                    e.define(t, i, {
                        configurable: !0,
                        get: function () {
                            return n[i];
                        }
                    });
                });
                for (var i in this)
                    !/^_/.test(i) && this[i] && (t[i] = this[i]);
            },
            setup: function (t) {
                return (paper = this), (this.project = new v(t)), this;
            },
            activate: function () {
                paper = this;
            },
            clear: function () {
                for (var t = this.projects.length - 1; t >= 0; t--)
                    this.projects[t].remove();
                for (var t = this.tools.length - 1; t >= 0; t--)
                    this.tools[t].remove();
                for (var t = this.palettes.length - 1; t >= 0; t--)
                    this.palettes[t].remove();
            },
            remove: function () {
                this.clear(), delete r._scopes[this._id];
            },
            statics: new (function () {
                function t(t) {
                    return (
                        (t += "Attribute"),
                        function (e, n) {
                            return e[t](n) || e[t]("data-paper-" + n);
                        }
                    );
                }
                return {
                    _scopes: {},
                    _id: 0,
                    get: function (t) {
                        return this._scopes[t] || null;
                    },
                    getAttribute: t("get"),
                    hasAttribute: t("has")
                };
            })()
        }),
        s = e.extend(n, {
            initialize: function (t) {
                (this._scope = paper),
                    (this._index = this._scope[this._list].push(this) - 1),
                    (t || !this._scope[this._reference]) && this.activate();
            },
            activate: function () {
                if (!this._scope) return !1;
                var t = this._scope[this._reference];
                return (
                    t && t !== this && t.fire("deactivate"),
                    (this._scope[this._reference] = this),
                    this.fire("activate", t),
                    !0
                );
            },
            isActive: function () {
                return this._scope[this._reference] === this;
            },
            remove: function () {
                return null == this._index
                    ? !1
                    : (e.splice(this._scope[this._list], null, this._index, 1),
                      this._scope[this._reference] == this &&
                          (this._scope[this._reference] = null),
                      (this._scope = null),
                      !0);
            }
        }),
        a = e.extend({
            initialize: function (t) {
                (this.precision = t || 5),
                    (this.multiplier = Math.pow(10, this.precision));
            },
            number: function (t) {
                return Math.round(t * this.multiplier) / this.multiplier;
            },
            pair: function (t, e, n) {
                return this.number(t) + (n || ",") + this.number(e);
            },
            point: function (t, e) {
                return this.number(t.x) + (e || ",") + this.number(t.y);
            },
            size: function (t, e) {
                return (
                    this.number(t.width) + (e || ",") + this.number(t.height)
                );
            },
            rectangle: function (t, e) {
                return this.point(t, e) + (e || ",") + this.size(t, e);
            }
        });
    a.instance = new a();
    var o = new (function () {
            function e(e, n, i) {
                var r = n === t,
                    s = n - c,
                    a = i + c,
                    o = 0;
                return function (t) {
                    return (
                        (r || (t > s && a > t)) &&
                            (e[o++] = n > t ? n : t > i ? i : t),
                        o
                    );
                };
            }
            var n = [
                    [0.5773502691896257],
                    [0, 0.7745966692414834],
                    [0.33998104358485626, 0.8611363115940526],
                    [0, 0.5384693101056831, 0.906179845938664],
                    [0.2386191860831969, 0.6612093864662645, 0.932469514203152],
                    [
                        0, 0.4058451513773972, 0.7415311855993945,
                        0.9491079123427585
                    ],
                    [
                        0.1834346424956498, 0.525532409916329,
                        0.7966664774136267, 0.9602898564975363
                    ],
                    [
                        0, 0.3242534234038089, 0.6133714327005904,
                        0.8360311073266358, 0.9681602395076261
                    ],
                    [
                        0.14887433898163122, 0.4333953941292472,
                        0.6794095682990244, 0.8650633666889845,
                        0.9739065285171717
                    ],
                    [
                        0, 0.26954315595234496, 0.5190961292068118,
                        0.7301520055740494, 0.8870625997680953,
                        0.978228658146057
                    ],
                    [
                        0.1252334085114689, 0.3678314989981802,
                        0.5873179542866175, 0.7699026741943047,
                        0.9041172563704749, 0.9815606342467192
                    ],
                    [
                        0, 0.2304583159551348, 0.44849275103644687,
                        0.6423493394403402, 0.8015780907333099,
                        0.9175983992229779, 0.9841830547185881
                    ],
                    [
                        0.10805494870734367, 0.31911236892788974,
                        0.5152486363581541, 0.6872929048116855,
                        0.827201315069765, 0.9284348836635735,
                        0.9862838086968123
                    ],
                    [
                        0, 0.20119409399743451, 0.3941513470775634,
                        0.5709721726085388, 0.7244177313601701,
                        0.8482065834104272, 0.937273392400706,
                        0.9879925180204854
                    ],
                    [
                        0.09501250983763744, 0.2816035507792589,
                        0.45801677765722737, 0.6178762444026438,
                        0.755404408355003, 0.8656312023878318,
                        0.9445750230732326, 0.9894009349916499
                    ]
                ],
                i = [
                    [1],
                    [0.8888888888888888, 0.5555555555555556],
                    [0.6521451548625461, 0.34785484513745385],
                    [
                        0.5688888888888889, 0.47862867049936647,
                        0.23692688505618908
                    ],
                    [
                        0.46791393457269104, 0.3607615730481386,
                        0.17132449237917036
                    ],
                    [
                        0.4179591836734694, 0.3818300505051189,
                        0.27970539148927664, 0.1294849661688697
                    ],
                    [
                        0.362683783378362, 0.31370664587788727,
                        0.22238103445337448, 0.10122853629037626
                    ],
                    [
                        0.3302393550012598, 0.31234707704000286,
                        0.26061069640293544, 0.1806481606948574,
                        0.08127438836157441
                    ],
                    [
                        0.29552422471475287, 0.26926671930999635,
                        0.21908636251598204, 0.1494513491505806,
                        0.06667134430868814
                    ],
                    [
                        0.2729250867779006, 0.26280454451024665,
                        0.23319376459199048, 0.18629021092773426,
                        0.1255803694649046, 0.05566856711617366
                    ],
                    [
                        0.24914704581340277, 0.2334925365383548,
                        0.20316742672306592, 0.16007832854334622,
                        0.10693932599531843, 0.04717533638651183
                    ],
                    [
                        0.2325515532308739, 0.22628318026289723,
                        0.2078160475368885, 0.17814598076194574,
                        0.13887351021978725, 0.09212149983772845,
                        0.04048400476531588
                    ],
                    [
                        0.2152638534631578, 0.2051984637212956,
                        0.18553839747793782, 0.15720316715819355,
                        0.12151857068790319, 0.08015808715976021,
                        0.03511946033175186
                    ],
                    [
                        0.2025782419255613, 0.19843148532711158,
                        0.1861610000155622, 0.16626920581699392,
                        0.13957067792615432, 0.10715922046717194,
                        0.07036604748810812, 0.03075324199611727
                    ],
                    [
                        0.1894506104550685, 0.18260341504492358,
                        0.16915651939500254, 0.14959598881657674,
                        0.12462897125553388, 0.09515851168249279,
                        0.062253523938647894, 0.027152459411754096
                    ]
                ],
                r = Math.abs,
                s = Math.sqrt,
                a = Math.pow,
                h = Math.cos,
                u = Math.PI,
                l = 1e-5,
                c = 1e-11;
            return {
                TOLERANCE: l,
                EPSILON: c,
                KAPPA: (4 * (s(2) - 1)) / 3,
                isZero: function (t) {
                    return r(t) <= c;
                },
                integrate: function (t, e, r, s) {
                    for (
                        var a = n[s - 2],
                            o = i[s - 2],
                            h = 0.5 * (r - e),
                            u = h + e,
                            l = 0,
                            c = (s + 1) >> 1,
                            d = 1 & s ? o[l++] * t(u) : 0;
                        c > l;

                    ) {
                        var f = h * a[l];
                        d += o[l++] * (t(u + f) + t(u - f));
                    }
                    return h * d;
                },
                findRoot: function (t, e, n, i, s, a, o) {
                    for (var h = 0; a > h; h++) {
                        var u = t(n),
                            l = u / e(n),
                            c = n - l;
                        if (r(l) < o) return c;
                        u > 0
                            ? ((s = n), (n = i >= c ? 0.5 * (i + s) : c))
                            : ((i = n), (n = c >= s ? 0.5 * (i + s) : c));
                    }
                    return n;
                },
                solveQuadratic: function (t, n, i, a, o, h) {
                    var u = e(a, o, h);
                    if (r(t) < c)
                        return r(n) >= c ? u(-i / n) : r(i) < c ? -1 : 0;
                    var l = n / (2 * t),
                        d = i / t,
                        f = l * l;
                    if (d - c > f) return 0;
                    var _ = f > d ? s(f - d) : 0,
                        g = u(_ - l);
                    return _ > 0 && (g = u(-_ - l)), g;
                },
                solveCubic: function (t, n, i, l, d, f, _) {
                    if (r(t) < c) return o.solveQuadratic(n, i, l, d, f, _);
                    (n /= t), (i /= t), (l /= t);
                    var g = e(d, f, _),
                        p = n * n,
                        v = (p - 3 * i) / 9,
                        m = (2 * p * n - 9 * n * i + 27 * l) / 54,
                        y = v * v * v,
                        w = m * m - y;
                    if (((n /= 3), r(w) < c)) {
                        if (r(m) < c) return g(-n);
                        var x = s(v),
                            b = m > 0 ? 1 : -1;
                        return g(2 * -b * x - n), g(b * x - n);
                    }
                    if (0 > w) {
                        var x = s(v),
                            C = Math.acos(m / (x * x * x)) / 3,
                            S = -2 * x,
                            P = (2 * u) / 3;
                        return (
                            g(S * h(C) - n),
                            g(S * h(C + P) - n),
                            g(S * h(C - P) - n)
                        );
                    }
                    var k = (m > 0 ? -1 : 1) * a(r(m) + s(w), 1 / 3);
                    return g(k + v / k - n);
                }
            };
        })(),
        h = e.extend(
            {
                _class: "Point",
                _readIndex: !0,
                initialize: function (t, e) {
                    var n = typeof t;
                    if ("number" === n) {
                        var i = "number" == typeof e;
                        (this.x = t),
                            (this.y = i ? e : t),
                            this.__read && (this.__read = i ? 2 : 1);
                    } else
                        "undefined" === n || null === t
                            ? ((this.x = this.y = 0),
                              this.__read && (this.__read = null === t ? 1 : 0))
                            : (Array.isArray(t)
                                  ? ((this.x = t[0]),
                                    (this.y = t.length > 1 ? t[1] : t[0]))
                                  : null != t.x
                                  ? ((this.x = t.x), (this.y = t.y))
                                  : null != t.width
                                  ? ((this.x = t.width), (this.y = t.height))
                                  : null != t.angle
                                  ? ((this.x = t.length),
                                    (this.y = 0),
                                    this.setAngle(t.angle))
                                  : ((this.x = this.y = 0),
                                    this.__read && (this.__read = 0)),
                              this.__read && (this.__read = 1));
                },
                set: function (t, e) {
                    return (this.x = t), (this.y = e), this;
                },
                equals: function (t) {
                    return (
                        this === t ||
                        (t &&
                            ((this.x === t.x && this.y === t.y) ||
                                (Array.isArray(t) &&
                                    this.x === t[0] &&
                                    this.y === t[1]))) ||
                        !1
                    );
                },
                clone: function () {
                    return new h(this.x, this.y);
                },
                toString: function () {
                    var t = a.instance;
                    return (
                        "{ x: " +
                        t.number(this.x) +
                        ", y: " +
                        t.number(this.y) +
                        " }"
                    );
                },
                _serialize: function (t) {
                    var e = t.formatter;
                    return [e.number(this.x), e.number(this.y)];
                },
                getLength: function () {
                    return Math.sqrt(this.x * this.x + this.y * this.y);
                },
                setLength: function (t) {
                    if (this.isZero()) {
                        var e = this._angle || 0;
                        this.set(Math.cos(e) * t, Math.sin(e) * t);
                    } else {
                        var n = t / this.getLength();
                        o.isZero(n) && this.getAngle(),
                            this.set(this.x * n, this.y * n);
                    }
                },
                getAngle: function () {
                    return (
                        (180 * this.getAngleInRadians.apply(this, arguments)) /
                        Math.PI
                    );
                },
                setAngle: function (t) {
                    this.setAngleInRadians.call(this, (t * Math.PI) / 180);
                },
                getAngleInDegrees: "#getAngle",
                setAngleInDegrees: "#setAngle",
                getAngleInRadians: function () {
                    if (arguments.length) {
                        var t = h.read(arguments),
                            e = this.getLength() * t.getLength();
                        if (o.isZero(e)) return 0 / 0;
                        var n = this.dot(t) / e;
                        return Math.acos(-1 > n ? -1 : n > 1 ? 1 : n);
                    }
                    return this.isZero()
                        ? this._angle || 0
                        : (this._angle = Math.atan2(this.y, this.x));
                },
                setAngleInRadians: function (t) {
                    if (((this._angle = t), !this.isZero())) {
                        var e = this.getLength();
                        this.set(Math.cos(t) * e, Math.sin(t) * e);
                    }
                },
                getQuadrant: function () {
                    return this.x >= 0
                        ? this.y >= 0
                            ? 1
                            : 4
                        : this.y >= 0
                        ? 2
                        : 3;
                }
            },
            {
                beans: !1,
                getDirectedAngle: function () {
                    var t = h.read(arguments);
                    return (
                        (180 * Math.atan2(this.cross(t), this.dot(t))) / Math.PI
                    );
                },
                getDistance: function () {
                    var t = h.read(arguments),
                        n = t.x - this.x,
                        i = t.y - this.y,
                        r = n * n + i * i,
                        s = e.read(arguments);
                    return s ? r : Math.sqrt(r);
                },
                normalize: function (e) {
                    e === t && (e = 1);
                    var n = this.getLength(),
                        i = 0 !== n ? e / n : 0,
                        r = new h(this.x * i, this.y * i);
                    return i >= 0 && (r._angle = this._angle), r;
                },
                rotate: function (t, e) {
                    if (0 === t) return this.clone();
                    t = (t * Math.PI) / 180;
                    var n = e ? this.subtract(e) : this,
                        i = Math.sin(t),
                        r = Math.cos(t);
                    return (
                        (n = new h(n.x * r - n.y * i, n.x * i + n.y * r)),
                        e ? n.add(e) : n
                    );
                },
                transform: function (t) {
                    return t ? t._transformPoint(this) : this;
                },
                add: function () {
                    var t = h.read(arguments);
                    return new h(this.x + t.x, this.y + t.y);
                },
                subtract: function () {
                    var t = h.read(arguments);
                    return new h(this.x - t.x, this.y - t.y);
                },
                multiply: function () {
                    var t = h.read(arguments);
                    return new h(this.x * t.x, this.y * t.y);
                },
                divide: function () {
                    var t = h.read(arguments);
                    return new h(this.x / t.x, this.y / t.y);
                },
                modulo: function () {
                    var t = h.read(arguments);
                    return new h(this.x % t.x, this.y % t.y);
                },
                negate: function () {
                    return new h(-this.x, -this.y);
                },
                isInside: function (t) {
                    return t.contains(this);
                },
                isClose: function (t, e) {
                    return this.getDistance(t) < e;
                },
                isColinear: function (t) {
                    return Math.abs(this.cross(t)) < 1e-5;
                },
                isOrthogonal: function (t) {
                    return Math.abs(this.dot(t)) < 1e-5;
                },
                isZero: function () {
                    return o.isZero(this.x) && o.isZero(this.y);
                },
                isNaN: function () {
                    return isNaN(this.x) || isNaN(this.y);
                },
                dot: function () {
                    var t = h.read(arguments);
                    return this.x * t.x + this.y * t.y;
                },
                cross: function () {
                    var t = h.read(arguments);
                    return this.x * t.y - this.y * t.x;
                },
                project: function () {
                    var t = h.read(arguments);
                    if (t.isZero()) return new h(0, 0);
                    var e = this.dot(t) / t.dot(t);
                    return new h(t.x * e, t.y * e);
                },
                statics: {
                    min: function () {
                        var t = h.read(arguments),
                            e = h.read(arguments);
                        return new h(Math.min(t.x, e.x), Math.min(t.y, e.y));
                    },
                    max: function () {
                        var t = h.read(arguments),
                            e = h.read(arguments);
                        return new h(Math.max(t.x, e.x), Math.max(t.y, e.y));
                    },
                    random: function () {
                        return new h(Math.random(), Math.random());
                    }
                }
            },
            e.each(
                ["round", "ceil", "floor", "abs"],
                function (t) {
                    var e = Math[t];
                    this[t] = function () {
                        return new h(e(this.x), e(this.y));
                    };
                },
                {}
            )
        ),
        u = h.extend({
            initialize: function (t, e, n, i) {
                (this._x = t),
                    (this._y = e),
                    (this._owner = n),
                    (this._setter = i);
            },
            set: function (t, e, n) {
                return (
                    (this._x = t),
                    (this._y = e),
                    n || this._owner[this._setter](this),
                    this
                );
            },
            getX: function () {
                return this._x;
            },
            setX: function (t) {
                (this._x = t), this._owner[this._setter](this);
            },
            getY: function () {
                return this._y;
            },
            setY: function (t) {
                (this._y = t), this._owner[this._setter](this);
            }
        }),
        c = e.extend(
            {
                _class: "Size",
                _readIndex: !0,
                initialize: function (t, e) {
                    var n = typeof t;
                    if ("number" === n) {
                        var i = "number" == typeof e;
                        (this.width = t),
                            (this.height = i ? e : t),
                            this.__read && (this.__read = i ? 2 : 1);
                    } else
                        "undefined" === n || null === t
                            ? ((this.width = this.height = 0),
                              this.__read && (this.__read = null === t ? 1 : 0))
                            : (Array.isArray(t)
                                  ? ((this.width = t[0]),
                                    (this.height = t.length > 1 ? t[1] : t[0]))
                                  : null != t.width
                                  ? ((this.width = t.width),
                                    (this.height = t.height))
                                  : null != t.x
                                  ? ((this.width = t.x), (this.height = t.y))
                                  : ((this.width = this.height = 0),
                                    this.__read && (this.__read = 0)),
                              this.__read && (this.__read = 1));
                },
                set: function (t, e) {
                    return (this.width = t), (this.height = e), this;
                },
                equals: function (t) {
                    return (
                        t === this ||
                        (t &&
                            ((this.width === t.width &&
                                this.height === t.height) ||
                                (Array.isArray(t) &&
                                    this.width === t[0] &&
                                    this.height === t[1]))) ||
                        !1
                    );
                },
                clone: function () {
                    return new c(this.width, this.height);
                },
                toString: function () {
                    var t = a.instance;
                    return (
                        "{ width: " +
                        t.number(this.width) +
                        ", height: " +
                        t.number(this.height) +
                        " }"
                    );
                },
                _serialize: function (t) {
                    var e = t.formatter;
                    return [e.number(this.width), e.number(this.height)];
                },
                add: function () {
                    var t = c.read(arguments);
                    return new c(this.width + t.width, this.height + t.height);
                },
                subtract: function () {
                    var t = c.read(arguments);
                    return new c(this.width - t.width, this.height - t.height);
                },
                multiply: function () {
                    var t = c.read(arguments);
                    return new c(this.width * t.width, this.height * t.height);
                },
                divide: function () {
                    var t = c.read(arguments);
                    return new c(this.width / t.width, this.height / t.height);
                },
                modulo: function () {
                    var t = c.read(arguments);
                    return new c(this.width % t.width, this.height % t.height);
                },
                negate: function () {
                    return new c(-this.width, -this.height);
                },
                isZero: function () {
                    return o.isZero(this.width) && o.isZero(this.height);
                },
                isNaN: function () {
                    return isNaN(this.width) || isNaN(this.height);
                },
                statics: {
                    min: function (t, e) {
                        return new c(
                            Math.min(t.width, e.width),
                            Math.min(t.height, e.height)
                        );
                    },
                    max: function (t, e) {
                        return new c(
                            Math.max(t.width, e.width),
                            Math.max(t.height, e.height)
                        );
                    },
                    random: function () {
                        return new c(Math.random(), Math.random());
                    }
                }
            },
            e.each(
                ["round", "ceil", "floor", "abs"],
                function (t) {
                    var e = Math[t];
                    this[t] = function () {
                        return new c(e(this.width), e(this.height));
                    };
                },
                {}
            )
        ),
        d = c.extend({
            initialize: function (t, e, n, i) {
                (this._width = t),
                    (this._height = e),
                    (this._owner = n),
                    (this._setter = i);
            },
            set: function (t, e, n) {
                return (
                    (this._width = t),
                    (this._height = e),
                    n || this._owner[this._setter](this),
                    this
                );
            },
            getWidth: function () {
                return this._width;
            },
            setWidth: function (t) {
                (this._width = t), this._owner[this._setter](this);
            },
            getHeight: function () {
                return this._height;
            },
            setHeight: function (t) {
                (this._height = t), this._owner[this._setter](this);
            }
        }),
        f = e.extend(
            {
                _class: "Rectangle",
                _readIndex: !0,
                beans: !0,
                initialize: function (n, i, r, s) {
                    var a = typeof n,
                        o = 0;
                    if (
                        ("number" === a
                            ? ((this.x = n),
                              (this.y = i),
                              (this.width = r),
                              (this.height = s),
                              (o = 4))
                            : "undefined" === a || null === n
                            ? ((this.x = this.y = this.width = this.height = 0),
                              (o = null === n ? 1 : 0))
                            : 1 === arguments.length &&
                              (Array.isArray(n)
                                  ? ((this.x = n[0]),
                                    (this.y = n[1]),
                                    (this.width = n[2]),
                                    (this.height = n[3]),
                                    (o = 1))
                                  : n.x !== t || n.width !== t
                                  ? ((this.x = n.x || 0),
                                    (this.y = n.y || 0),
                                    (this.width = n.width || 0),
                                    (this.height = n.height || 0),
                                    (o = 1))
                                  : n.from === t &&
                                    n.to === t &&
                                    ((this.x =
                                        this.y =
                                        this.width =
                                        this.height =
                                            0),
                                    this._set(n),
                                    (o = 1))),
                        !o)
                    ) {
                        var u = h.readNamed(arguments, "from"),
                            l = e.peek(arguments);
                        if (
                            ((this.x = u.x),
                            (this.y = u.y),
                            (l && l.x !== t) || e.hasNamed(arguments, "to"))
                        ) {
                            var d = h.readNamed(arguments, "to");
                            (this.width = d.x - u.x),
                                (this.height = d.y - u.y),
                                this.width < 0 &&
                                    ((this.x = d.x),
                                    (this.width = -this.width)),
                                this.height < 0 &&
                                    ((this.y = d.y),
                                    (this.height = -this.height));
                        } else {
                            var f = c.read(arguments);
                            (this.width = f.width), (this.height = f.height);
                        }
                        o = arguments.__index;
                    }
                    this.__read && (this.__read = o);
                },
                set: function (t, e, n, i) {
                    return (
                        (this.x = t),
                        (this.y = e),
                        (this.width = n),
                        (this.height = i),
                        this
                    );
                },
                clone: function () {
                    return new f(this.x, this.y, this.width, this.height);
                },
                equals: function (t) {
                    var n = e.isPlainValue(t) ? f.read(arguments) : t;
                    return (
                        n === this ||
                        (n &&
                            this.x === n.x &&
                            this.y === n.y &&
                            this.width === n.width &&
                            this.height === n.height) ||
                        !1
                    );
                },
                toString: function () {
                    var t = a.instance;
                    return (
                        "{ x: " +
                        t.number(this.x) +
                        ", y: " +
                        t.number(this.y) +
                        ", width: " +
                        t.number(this.width) +
                        ", height: " +
                        t.number(this.height) +
                        " }"
                    );
                },
                _serialize: function (t) {
                    var e = t.formatter;
                    return [
                        e.number(this.x),
                        e.number(this.y),
                        e.number(this.width),
                        e.number(this.height)
                    ];
                },
                getPoint: function (t) {
                    var e = t ? h : u;
                    return new e(this.x, this.y, this, "setPoint");
                },
                setPoint: function () {
                    var t = h.read(arguments);
                    (this.x = t.x), (this.y = t.y);
                },
                getSize: function (t) {
                    var e = t ? c : d;
                    return new e(this.width, this.height, this, "setSize");
                },
                setSize: function () {
                    var t = c.read(arguments);
                    this._fixX &&
                        (this.x += (this.width - t.width) * this._fixX),
                        this._fixY &&
                            (this.y += (this.height - t.height) * this._fixY),
                        (this.width = t.width),
                        (this.height = t.height),
                        (this._fixW = 1),
                        (this._fixH = 1);
                },
                getLeft: function () {
                    return this.x;
                },
                setLeft: function (t) {
                    this._fixW || (this.width -= t - this.x),
                        (this.x = t),
                        (this._fixX = 0);
                },
                getTop: function () {
                    return this.y;
                },
                setTop: function (t) {
                    this._fixH || (this.height -= t - this.y),
                        (this.y = t),
                        (this._fixY = 0);
                },
                getRight: function () {
                    return this.x + this.width;
                },
                setRight: function (e) {
                    this._fixX !== t && 1 !== this._fixX && (this._fixW = 0),
                        this._fixW
                            ? (this.x = e - this.width)
                            : (this.width = e - this.x),
                        (this._fixX = 1);
                },
                getBottom: function () {
                    return this.y + this.height;
                },
                setBottom: function (e) {
                    this._fixY !== t && 1 !== this._fixY && (this._fixH = 0),
                        this._fixH
                            ? (this.y = e - this.height)
                            : (this.height = e - this.y),
                        (this._fixY = 1);
                },
                getCenterX: function () {
                    return this.x + 0.5 * this.width;
                },
                setCenterX: function (t) {
                    (this.x = t - 0.5 * this.width), (this._fixX = 0.5);
                },
                getCenterY: function () {
                    return this.y + 0.5 * this.height;
                },
                setCenterY: function (t) {
                    (this.y = t - 0.5 * this.height), (this._fixY = 0.5);
                },
                getCenter: function (t) {
                    var e = t ? h : u;
                    return new e(
                        this.getCenterX(),
                        this.getCenterY(),
                        this,
                        "setCenter"
                    );
                },
                setCenter: function () {
                    var t = h.read(arguments);
                    return this.setCenterX(t.x), this.setCenterY(t.y), this;
                },
                getArea: function () {
                    return this.width * this.height;
                },
                isEmpty: function () {
                    return 0 === this.width || 0 === this.height;
                },
                contains: function (e) {
                    return (e && e.width !== t) ||
                        4 == (Array.isArray(e) ? e : arguments).length
                        ? this._containsRectangle(f.read(arguments))
                        : this._containsPoint(h.read(arguments));
                },
                _containsPoint: function (t) {
                    var e = t.x,
                        n = t.y;
                    return (
                        e >= this.x &&
                        n >= this.y &&
                        e <= this.x + this.width &&
                        n <= this.y + this.height
                    );
                },
                _containsRectangle: function (t) {
                    var e = t.x,
                        n = t.y;
                    return (
                        e >= this.x &&
                        n >= this.y &&
                        e + t.width <= this.x + this.width &&
                        n + t.height <= this.y + this.height
                    );
                },
                intersects: function () {
                    var t = f.read(arguments);
                    return (
                        t.x + t.width > this.x &&
                        t.y + t.height > this.y &&
                        t.x < this.x + this.width &&
                        t.y < this.y + this.height
                    );
                },
                touches: function () {
                    var t = f.read(arguments);
                    return (
                        t.x + t.width >= this.x &&
                        t.y + t.height >= this.y &&
                        t.x <= this.x + this.width &&
                        t.y <= this.y + this.height
                    );
                },
                intersect: function () {
                    var t = f.read(arguments),
                        e = Math.max(this.x, t.x),
                        n = Math.max(this.y, t.y),
                        i = Math.min(this.x + this.width, t.x + t.width),
                        r = Math.min(this.y + this.height, t.y + t.height);
                    return new f(e, n, i - e, r - n);
                },
                unite: function () {
                    var t = f.read(arguments),
                        e = Math.min(this.x, t.x),
                        n = Math.min(this.y, t.y),
                        i = Math.max(this.x + this.width, t.x + t.width),
                        r = Math.max(this.y + this.height, t.y + t.height);
                    return new f(e, n, i - e, r - n);
                },
                include: function () {
                    var t = h.read(arguments),
                        e = Math.min(this.x, t.x),
                        n = Math.min(this.y, t.y),
                        i = Math.max(this.x + this.width, t.x),
                        r = Math.max(this.y + this.height, t.y);
                    return new f(e, n, i - e, r - n);
                },
                expand: function () {
                    var t = c.read(arguments),
                        e = t.width,
                        n = t.height;
                    return new f(
                        this.x - e / 2,
                        this.y - n / 2,
                        this.width + e,
                        this.height + n
                    );
                },
                scale: function (e, n) {
                    return this.expand(
                        this.width * e - this.width,
                        this.height * (n === t ? e : n) - this.height
                    );
                }
            },
            e.each(
                [
                    ["Top", "Left"],
                    ["Top", "Right"],
                    ["Bottom", "Left"],
                    ["Bottom", "Right"],
                    ["Left", "Center"],
                    ["Top", "Center"],
                    ["Right", "Center"],
                    ["Bottom", "Center"]
                ],
                function (t, e) {
                    var n = t.join(""),
                        i = /^[RL]/.test(n);
                    e >= 4 && (t[1] += i ? "Y" : "X");
                    var r = t[i ? 0 : 1],
                        s = t[i ? 1 : 0],
                        a = "get" + r,
                        o = "get" + s,
                        l = "set" + r,
                        c = "set" + s,
                        d = "get" + n,
                        f = "set" + n;
                    (this[d] = function (t) {
                        var e = t ? h : u;
                        return new e(this[a](), this[o](), this, f);
                    }),
                        (this[f] = function () {
                            var t = h.read(arguments);
                            this[l](t.x), this[c](t.y);
                        });
                },
                {
                    beans: !0
                }
            )
        ),
        _ = f.extend(
            {
                initialize: function (t, e, n, i, r, s) {
                    this.set(t, e, n, i, !0),
                        (this._owner = r),
                        (this._setter = s);
                },
                set: function (t, e, n, i, r) {
                    return (
                        (this._x = t),
                        (this._y = e),
                        (this._width = n),
                        (this._height = i),
                        r || this._owner[this._setter](this),
                        this
                    );
                }
            },
            new (function () {
                var t = f.prototype;
                return e.each(
                    ["x", "y", "width", "height"],
                    function (t) {
                        var n = e.capitalize(t),
                            i = "_" + t;
                        (this["get" + n] = function () {
                            return this[i];
                        }),
                            (this["set" + n] = function (t) {
                                (this[i] = t),
                                    this._dontNotify ||
                                        this._owner[this._setter](this);
                            });
                    },
                    e.each(
                        [
                            "Point",
                            "Size",
                            "Center",
                            "Left",
                            "Top",
                            "Right",
                            "Bottom",
                            "CenterX",
                            "CenterY",
                            "TopLeft",
                            "TopRight",
                            "BottomLeft",
                            "BottomRight",
                            "LeftCenter",
                            "TopCenter",
                            "RightCenter",
                            "BottomCenter"
                        ],
                        function (e) {
                            var n = "set" + e;
                            this[n] = function () {
                                (this._dontNotify = !0),
                                    t[n].apply(this, arguments),
                                    (this._dontNotify = !1),
                                    this._owner[this._setter](this);
                            };
                        },
                        {
                            isSelected: function () {
                                return this._owner._boundsSelected;
                            },
                            setSelected: function (t) {
                                var e = this._owner;
                                e.setSelected &&
                                    ((e._boundsSelected = t),
                                    e.setSelected(
                                        t || e._selectedSegmentState > 0
                                    ));
                            }
                        }
                    )
                );
            })()
        ),
        g = e.extend(
            {
                _class: "Matrix",
                initialize: function ae(t) {
                    var e = arguments.length,
                        n = !0;
                    if (
                        (6 === e
                            ? this.set.apply(this, arguments)
                            : 1 === e
                            ? t instanceof ae
                                ? this.set(t._a, t._c, t._b, t._d, t._tx, t._ty)
                                : Array.isArray(t)
                                ? this.set.apply(this, t)
                                : (n = !1)
                            : 0 === e
                            ? this.reset()
                            : (n = !1),
                        !n)
                    )
                        throw Error("Unsupported matrix parameters");
                },
                set: function (t, e, n, i, r, s, a) {
                    return (
                        (this._a = t),
                        (this._c = e),
                        (this._b = n),
                        (this._d = i),
                        (this._tx = r),
                        (this._ty = s),
                        a || this._changed(),
                        this
                    );
                },
                _serialize: function (t) {
                    return e.serialize(this.getValues(), t);
                },
                _changed: function () {
                    var t = this._owner;
                    t &&
                        (t._applyMatrix
                            ? t.transform(null, !0)
                            : t._changed(9));
                },
                clone: function () {
                    return new g(
                        this._a,
                        this._c,
                        this._b,
                        this._d,
                        this._tx,
                        this._ty
                    );
                },
                equals: function (t) {
                    return (
                        t === this ||
                        (t &&
                            this._a === t._a &&
                            this._b === t._b &&
                            this._c === t._c &&
                            this._d === t._d &&
                            this._tx === t._tx &&
                            this._ty === t._ty) ||
                        !1
                    );
                },
                toString: function () {
                    var t = a.instance;
                    return (
                        "[[" +
                        [
                            t.number(this._a),
                            t.number(this._b),
                            t.number(this._tx)
                        ].join(", ") +
                        "], [" +
                        [
                            t.number(this._c),
                            t.number(this._d),
                            t.number(this._ty)
                        ].join(", ") +
                        "]]"
                    );
                },
                reset: function (t) {
                    return (
                        (this._a = this._d = 1),
                        (this._c = this._b = this._tx = this._ty = 0),
                        t || this._changed(),
                        this
                    );
                },
                apply: function () {
                    var t = this._owner;
                    return t ? (t.transform(null, !0), this.isIdentity()) : !1;
                },
                translate: function () {
                    var t = h.read(arguments),
                        e = t.x,
                        n = t.y;
                    return (
                        (this._tx += e * this._a + n * this._b),
                        (this._ty += e * this._c + n * this._d),
                        this._changed(),
                        this
                    );
                },
                scale: function () {
                    var t = h.read(arguments),
                        e = h.read(arguments, 0, {
                            readNull: !0
                        });
                    return (
                        e && this.translate(e),
                        (this._a *= t.x),
                        (this._c *= t.x),
                        (this._b *= t.y),
                        (this._d *= t.y),
                        e && this.translate(e.negate()),
                        this._changed(),
                        this
                    );
                },
                rotate: function (t) {
                    t *= Math.PI / 180;
                    var e = h.read(arguments, 1),
                        n = e.x,
                        i = e.y,
                        r = Math.cos(t),
                        s = Math.sin(t),
                        a = n - n * r + i * s,
                        o = i - n * s - i * r,
                        u = this._a,
                        l = this._b,
                        c = this._c,
                        d = this._d;
                    return (
                        (this._a = r * u + s * l),
                        (this._b = -s * u + r * l),
                        (this._c = r * c + s * d),
                        (this._d = -s * c + r * d),
                        (this._tx += a * u + o * l),
                        (this._ty += a * c + o * d),
                        this._changed(),
                        this
                    );
                },
                shear: function () {
                    var t = h.read(arguments),
                        e = h.read(arguments, 0, {
                            readNull: !0
                        });
                    e && this.translate(e);
                    var n = this._a,
                        i = this._c;
                    return (
                        (this._a += t.y * this._b),
                        (this._c += t.y * this._d),
                        (this._b += t.x * n),
                        (this._d += t.x * i),
                        e && this.translate(e.negate()),
                        this._changed(),
                        this
                    );
                },
                skew: function () {
                    var t = h.read(arguments),
                        e = h.read(arguments, 0, {
                            readNull: !0
                        }),
                        n = Math.PI / 180,
                        i = new h(Math.tan(t.x * n), Math.tan(t.y * n));
                    return this.shear(i, e);
                },
                concatenate: function (t) {
                    var e = this._a,
                        n = this._b,
                        i = this._c,
                        r = this._d,
                        s = t._a,
                        a = t._b,
                        o = t._c,
                        h = t._d,
                        u = t._tx,
                        l = t._ty;
                    return (
                        (this._a = s * e + o * n),
                        (this._b = a * e + h * n),
                        (this._c = s * i + o * r),
                        (this._d = a * i + h * r),
                        (this._tx += u * e + l * n),
                        (this._ty += u * i + l * r),
                        this._changed(),
                        this
                    );
                },
                preConcatenate: function (t) {
                    var e = this._a,
                        n = this._b,
                        i = this._c,
                        r = this._d,
                        s = this._tx,
                        a = this._ty,
                        o = t._a,
                        h = t._b,
                        u = t._c,
                        l = t._d,
                        c = t._tx,
                        d = t._ty;
                    return (
                        (this._a = o * e + h * i),
                        (this._b = o * n + h * r),
                        (this._c = u * e + l * i),
                        (this._d = u * n + l * r),
                        (this._tx = o * s + h * a + c),
                        (this._ty = u * s + l * a + d),
                        this._changed(),
                        this
                    );
                },
                chain: function (t) {
                    var e = this._a,
                        n = this._b,
                        i = this._c,
                        r = this._d,
                        s = this._tx,
                        a = this._ty,
                        o = t._a,
                        h = t._b,
                        u = t._c,
                        l = t._d,
                        c = t._tx,
                        d = t._ty;
                    return new g(
                        o * e + u * n,
                        o * i + u * r,
                        h * e + l * n,
                        h * i + l * r,
                        s + c * e + d * n,
                        a + c * i + d * r
                    );
                },
                isIdentity: function () {
                    return (
                        1 === this._a &&
                        0 === this._c &&
                        0 === this._b &&
                        1 === this._d &&
                        0 === this._tx &&
                        0 === this._ty
                    );
                },
                orNullIfIdentity: function () {
                    return this.isIdentity() ? null : this;
                },
                isInvertible: function () {
                    return !!this._getDeterminant();
                },
                isSingular: function () {
                    return !this._getDeterminant();
                },
                transform: function (t, e, n) {
                    return arguments.length < 3
                        ? this._transformPoint(h.read(arguments))
                        : this._transformCoordinates(t, e, n);
                },
                _transformPoint: function (t, e, n) {
                    var i = t.x,
                        r = t.y;
                    return (
                        e || (e = new h()),
                        e.set(
                            i * this._a + r * this._b + this._tx,
                            i * this._c + r * this._d + this._ty,
                            n
                        )
                    );
                },
                _transformCoordinates: function (t, e, n) {
                    for (var i = 0, r = 0, s = 2 * n; s > i; ) {
                        var a = t[i++],
                            o = t[i++];
                        (e[r++] = a * this._a + o * this._b + this._tx),
                            (e[r++] = a * this._c + o * this._d + this._ty);
                    }
                    return e;
                },
                _transformCorners: function (t) {
                    var e = t.x,
                        n = t.y,
                        i = e + t.width,
                        r = n + t.height,
                        s = [e, n, i, n, i, r, e, r];
                    return this._transformCoordinates(s, s, 4);
                },
                _transformBounds: function (t, e, n) {
                    for (
                        var i = this._transformCorners(t),
                            r = i.slice(0, 2),
                            s = i.slice(),
                            a = 2;
                        8 > a;
                        a++
                    ) {
                        var o = i[a],
                            h = 1 & a;
                        o < r[h] ? (r[h] = o) : o > s[h] && (s[h] = o);
                    }
                    return (
                        e || (e = new f()),
                        e.set(r[0], r[1], s[0] - r[0], s[1] - r[1], n)
                    );
                },
                inverseTransform: function () {
                    return this._inverseTransform(h.read(arguments));
                },
                _getDeterminant: function () {
                    var t = this._a * this._d - this._b * this._c;
                    return isFinite(t) &&
                        !o.isZero(t) &&
                        isFinite(this._tx) &&
                        isFinite(this._ty)
                        ? t
                        : null;
                },
                _inverseTransform: function (t, e, n) {
                    var i = this._getDeterminant();
                    if (!i) return null;
                    var r = t.x - this._tx,
                        s = t.y - this._ty;
                    return (
                        e || (e = new h()),
                        e.set(
                            (r * this._d - s * this._b) / i,
                            (s * this._a - r * this._c) / i,
                            n
                        )
                    );
                },
                decompose: function () {
                    var t = this._a,
                        e = this._b,
                        n = this._c,
                        i = this._d;
                    if (o.isZero(t * i - e * n)) return null;
                    var r = Math.sqrt(t * t + e * e);
                    (t /= r), (e /= r);
                    var s = t * n + e * i;
                    (n -= t * s), (i -= e * s);
                    var a = Math.sqrt(n * n + i * i);
                    return (
                        (n /= a),
                        (i /= a),
                        (s /= a),
                        e * n > t * i &&
                            ((t = -t), (e = -e), (s = -s), (r = -r)),
                        {
                            scaling: new h(r, a),
                            rotation: (180 * -Math.atan2(e, t)) / Math.PI,
                            shearing: s
                        }
                    );
                },
                getValues: function () {
                    return [
                        this._a,
                        this._c,
                        this._b,
                        this._d,
                        this._tx,
                        this._ty
                    ];
                },
                getTranslation: function () {
                    return new h(this._tx, this._ty);
                },
                getScaling: function () {
                    return (this.decompose() || {}).scaling;
                },
                getRotation: function () {
                    return (this.decompose() || {}).rotation;
                },
                inverted: function () {
                    var t = this._getDeterminant();
                    return (
                        t &&
                        new g(
                            this._d / t,
                            -this._c / t,
                            -this._b / t,
                            this._a / t,
                            (this._b * this._ty - this._d * this._tx) / t,
                            (this._c * this._tx - this._a * this._ty) / t
                        )
                    );
                },
                shiftless: function () {
                    return new g(this._a, this._c, this._b, this._d, 0, 0);
                },
                applyToContext: function (t) {
                    t.transform(
                        this._a,
                        this._c,
                        this._b,
                        this._d,
                        this._tx,
                        this._ty
                    );
                }
            },
            e.each(
                ["a", "c", "b", "d", "tx", "ty"],
                function (t) {
                    var n = e.capitalize(t),
                        i = "_" + t;
                    (this["get" + n] = function () {
                        return this[i];
                    }),
                        (this["set" + n] = function (t) {
                            (this[i] = t), this._changed();
                        });
                },
                {}
            )
        ),
        p = e.extend({
            _class: "Line",
            initialize: function (t, e, n, i, r) {
                var s = !1;
                arguments.length >= 4
                    ? ((this._px = t),
                      (this._py = e),
                      (this._vx = n),
                      (this._vy = i),
                      (s = r))
                    : ((this._px = t.x),
                      (this._py = t.y),
                      (this._vx = e.x),
                      (this._vy = e.y),
                      (s = n)),
                    s || ((this._vx -= this._px), (this._vy -= this._py));
            },
            getPoint: function () {
                return new h(this._px, this._py);
            },
            getVector: function () {
                return new h(this._vx, this._vy);
            },
            getLength: function () {
                return this.getVector().getLength();
            },
            intersect: function (t, e) {
                return p.intersect(
                    this._px,
                    this._py,
                    this._vx,
                    this._vy,
                    t._px,
                    t._py,
                    t._vx,
                    t._vy,
                    !0,
                    e
                );
            },
            getSide: function (t) {
                return p.getSide(
                    this._px,
                    this._py,
                    this._vx,
                    this._vy,
                    t.x,
                    t.y,
                    !0
                );
            },
            getDistance: function (t) {
                return Math.abs(
                    p.getSignedDistance(
                        this._px,
                        this._py,
                        this._vx,
                        this._vy,
                        t.x,
                        t.y,
                        !0
                    )
                );
            },
            statics: {
                intersect: function (t, e, n, i, r, s, a, u, l, c) {
                    l || ((n -= t), (i -= e), (a -= r), (u -= s));
                    var d = u * n - a * i;
                    if (!o.isZero(d)) {
                        var f = t - r,
                            _ = e - s,
                            g = (a * _ - u * f) / d,
                            p = (n * _ - i * f) / d;
                        if (
                            (c || (g >= 0 && 1 >= g)) &&
                            (c || (p >= 0 && 1 >= p))
                        )
                            return new h(t + g * n, e + g * i);
                    }
                },
                getSide: function (t, e, n, i, r, s, a) {
                    a || ((n -= t), (i -= e));
                    var o = r - t,
                        h = s - e,
                        u = o * i - h * n;
                    return (
                        0 === u &&
                            ((u = o * n + h * i),
                            u > 0 &&
                                ((o -= n),
                                (h -= i),
                                (u = o * n + h * i),
                                0 > u && (u = 0))),
                        0 > u ? -1 : u > 0 ? 1 : 0
                    );
                },
                getSignedDistance: function (t, e, n, i, r, s, a) {
                    a || ((n -= t), (i -= e));
                    var o = i / n,
                        h = e - o * t;
                    return (s - o * r - h) / Math.sqrt(o * o + 1);
                }
            }
        }),
        v = s.extend({
            _class: "Project",
            _list: "projects",
            _reference: "project",
            initialize: function (t) {
                s.call(this, !0),
                    (this.layers = []),
                    (this.symbols = []),
                    (this._currentStyle = new F(null, null, this)),
                    (this.activeLayer = new x()),
                    (this._view = Z.create(this, t || Q.getCanvas(1, 1))),
                    (this._selectedItems = {}),
                    (this._selectedItemCount = 0),
                    (this._updateVersion = 0);
            },
            _serialize: function (t, n) {
                return e.serialize(this.layers, t, !0, n);
            },
            clear: function () {
                for (var t = this.layers.length - 1; t >= 0; t--)
                    this.layers[t].remove();
                this.symbols = [];
            },
            isEmpty: function () {
                return (
                    this.layers.length <= 1 &&
                    (!this.activeLayer || this.activeLayer.isEmpty())
                );
            },
            remove: function oe() {
                return oe.base.call(this)
                    ? (this._view && this._view.remove(), !0)
                    : !1;
            },
            getView: function () {
                return this._view;
            },
            getCurrentStyle: function () {
                return this._currentStyle;
            },
            setCurrentStyle: function (t) {
                this._currentStyle.initialize(t);
            },
            getIndex: function () {
                return this._index;
            },
            addChild: function (t) {
                return (
                    t instanceof x
                        ? (e.splice(this.layers, [t]),
                          this.activeLayer || (this.activeLayer = t))
                        : t instanceof y
                        ? (
                              this.activeLayer ||
                              this.addChild(new x(y.NO_INSERT))
                          ).addChild(t)
                        : (t = null),
                    t
                );
            },
            getSelectedItems: function () {
                var t = [];
                for (var e in this._selectedItems) {
                    var n = this._selectedItems[e];
                    n.isInserted() && t.push(n);
                }
                return t;
            },
            getOptions: function () {
                return this._scope.settings;
            },
            _updateSelection: function (t) {
                var e = t._id,
                    n = this._selectedItems;
                t._selected
                    ? n[e] !== t && (this._selectedItemCount++, (n[e] = t))
                    : n[e] === t && (this._selectedItemCount--, delete n[e]);
            },
            selectAll: function () {
                for (var t = this.layers, e = 0, n = t.length; n > e; e++)
                    t[e].setFullySelected(!0);
            },
            deselectAll: function () {
                var t = this._selectedItems;
                for (var e in t) t[e].setFullySelected(!1);
            },
            hitTest: function () {
                for (
                    var t = h.read(arguments),
                        n = P.getOptions(e.read(arguments)),
                        i = this.layers.length - 1;
                    i >= 0;
                    i--
                ) {
                    var r = this.layers[i]._hitTest(t, n);
                    if (r) return r;
                }
                return null;
            },
            getItems: function (t) {
                return y._getItems(this.layers, t, !0);
            },
            getItem: function (t) {
                return y._getItems(this.layers, t, !1);
            },
            importJSON: function (t) {
                this.activate();
                var n = this.activeLayer;
                return e.importJSON(t, n && n.isEmpty() && n);
            },
            draw: function (t, n, i) {
                this._updateVersion++, t.save(), n.applyToContext(t);
                for (
                    var r = new e({
                            offset: new h(0, 0),
                            pixelRatio: i,
                            viewMatrix: n.isIdentity() ? null : n,
                            matrices: [new g()],
                            updateMatrix: !0
                        }),
                        s = 0,
                        a = this.layers,
                        o = a.length;
                    o > s;
                    s++
                )
                    a[s].draw(t, r);
                if ((t.restore(), this._selectedItemCount > 0)) {
                    t.save(), (t.strokeWidth = 1);
                    var u = this._selectedItems,
                        l = this._scope.settings.handleSize,
                        c = this._updateVersion;
                    for (var d in u) u[d]._drawSelection(t, n, l, u, c);
                    t.restore();
                }
            }
        }),
        m = e.extend({
            _class: "Symbol",
            initialize: function he(t, e) {
                (this._id = he._id = (he._id || 0) + 1),
                    (this.project = paper.project),
                    this.project.symbols.push(this),
                    t && this.setDefinition(t, e);
            },
            _serialize: function (t, n) {
                return n.add(this, function () {
                    return e.serialize(
                        [this._class, this._definition],
                        t,
                        !1,
                        n
                    );
                });
            },
            _changed: function (t) {
                8 & t && y._clearBoundsCache(this),
                    1 & t && (this.project._needsUpdate = !0);
            },
            getDefinition: function () {
                return this._definition;
            },
            setDefinition: function (t, e) {
                t._parentSymbol && (t = t.clone()),
                    this._definition && (this._definition._parentSymbol = null),
                    (this._definition = t),
                    t.remove(),
                    t.setSelected(!1),
                    e || t.setPosition(new h()),
                    (t._parentSymbol = this),
                    this._changed(9);
            },
            place: function (t) {
                return new S(this, t);
            },
            clone: function () {
                return new m(this._definition.clone(!1));
            }
        }),
        y = e.extend(
            n,
            {
                statics: {
                    extend: function ue(t) {
                        return (
                            t._serializeFields &&
                                (t._serializeFields = new e(
                                    this.prototype._serializeFields,
                                    t._serializeFields
                                )),
                            ue.base.apply(this, arguments)
                        );
                    },
                    NO_INSERT: {
                        insert: !1
                    }
                },
                _class: "Item",
                _applyMatrix: !0,
                _canApplyMatrix: !0,
                _boundsSelected: !1,
                _selectChildren: !1,
                _serializeFields: {
                    name: null,
                    applyMatrix: null,
                    matrix: new g(),
                    pivot: null,
                    locked: !1,
                    visible: !0,
                    blendMode: "normal",
                    opacity: 1,
                    guide: !1,
                    selected: !1,
                    clipMask: !1,
                    data: {}
                },
                initialize: function () {},
                _initialize: function (t, n) {
                    var i = t && e.isPlainObject(t),
                        r = i && t.internal === !0,
                        s = (this._matrix = new g()),
                        a = paper.project;
                    return (
                        r || (this._id = y._id = (y._id || 0) + 1),
                        (this._applyMatrix =
                            this._canApplyMatrix && paper.settings.applyMatrix),
                        n && s.translate(n),
                        (s._owner = this),
                        (this._style = new F(a._currentStyle, this, a)),
                        this._project ||
                            (r || (i && t.insert === !1)
                                ? this._setProject(a)
                                : i && t.parent
                                ? this.setParent(t.parent)
                                : (a.activeLayer || new x()).addChild(this)),
                        i &&
                            t !== y.NO_INSERT &&
                            this._set(
                                t,
                                {
                                    insert: !0,
                                    parent: !0
                                },
                                !0
                            ),
                        i
                    );
                },
                _events: new (function () {
                    var t = {
                            mousedown: {
                                mousedown: 1,
                                mousedrag: 1,
                                click: 1,
                                doubleclick: 1
                            },
                            mouseup: {
                                mouseup: 1,
                                mousedrag: 1,
                                click: 1,
                                doubleclick: 1
                            },
                            mousemove: {
                                mousedrag: 1,
                                mousemove: 1,
                                mouseenter: 1,
                                mouseleave: 1
                            }
                        },
                        n = {
                            install: function (e) {
                                var n = this.getView()._eventCounters;
                                if (n)
                                    for (var i in t)
                                        n[i] = (n[i] || 0) + (t[i][e] || 0);
                            },
                            uninstall: function (e) {
                                var n = this.getView()._eventCounters;
                                if (n) for (var i in t) n[i] -= t[i][e] || 0;
                            }
                        };
                    return e.each(
                        [
                            "onMouseDown",
                            "onMouseUp",
                            "onMouseDrag",
                            "onClick",
                            "onDoubleClick",
                            "onMouseMove",
                            "onMouseEnter",
                            "onMouseLeave"
                        ],
                        function (t) {
                            this[t] = n;
                        },
                        {
                            onFrame: {
                                install: function () {
                                    this._animateItem(!0);
                                },
                                uninstall: function () {
                                    this._animateItem(!1);
                                }
                            },
                            onLoad: {}
                        }
                    );
                })(),
                _animateItem: function (t) {
                    this.getView()._animateItem(this, t);
                },
                _serialize: function (t, n) {
                    function i(i) {
                        for (var a in i) {
                            var o = s[a];
                            e.equals(
                                o,
                                "leading" === a ? 1.2 * i.fontSize : i[a]
                            ) || (r[a] = e.serialize(o, t, "data" !== a, n));
                        }
                    }
                    var r = {},
                        s = this;
                    return (
                        i(this._serializeFields),
                        this instanceof w || i(this._style._defaults),
                        [this._class, r]
                    );
                },
                _changed: function (e) {
                    var n = this._parentSymbol,
                        i = this._parent || n,
                        r = this._project;
                    if (
                        (8 & e &&
                            (this._bounds =
                                this._position =
                                this._decomposed =
                                this._globalMatrix =
                                this._currentPath =
                                    t),
                        i && 40 & e && y._clearBoundsCache(i),
                        2 & e && y._clearBoundsCache(this),
                        r && (1 & e && (r._needsUpdate = !0), r._changes))
                    ) {
                        var s = r._changesById[this._id];
                        s
                            ? (s.flags |= e)
                            : ((s = {
                                  item: this,
                                  flags: e
                              }),
                              (r._changesById[this._id] = s),
                              r._changes.push(s));
                    }
                    n && n._changed(e);
                },
                set: function (t) {
                    return t && this._set(t), this;
                },
                getId: function () {
                    return this._id;
                },
                getClassName: function () {
                    return this._class;
                },
                getName: function () {
                    return this._name;
                },
                setName: function (e, n) {
                    if ((this._name && this._removeNamed(), e === +e + ""))
                        throw Error(
                            "Names consisting only of numbers are not supported."
                        );
                    var i = this._parent;
                    if (e && i) {
                        for (
                            var r = i._children,
                                s = i._namedChildren,
                                a = e,
                                o = 1;
                            n && r[e];

                        )
                            e = a + " " + o++;
                        (s[e] = s[e] || []).push(this), (r[e] = this);
                    }
                    (this._name = e || t), this._changed(128);
                },
                getStyle: function () {
                    return this._style;
                },
                setStyle: function (t) {
                    this.getStyle().set(t);
                }
            },
            e.each(
                ["locked", "visible", "blendMode", "opacity", "guide"],
                function (t) {
                    var n = e.capitalize(t),
                        t = "_" + t;
                    (this["get" + n] = function () {
                        return this[t];
                    }),
                        (this["set" + n] = function (e) {
                            e != this[t] &&
                                ((this[t] = e),
                                this._changed("_locked" === t ? 128 : 129));
                        });
                },
                {}
            ),
            {
                beans: !0,
                _locked: !1,
                _visible: !0,
                _blendMode: "normal",
                _opacity: 1,
                _guide: !1,
                isSelected: function () {
                    if (this._selectChildren)
                        for (
                            var t = this._children, e = 0, n = t.length;
                            n > e;
                            e++
                        )
                            if (t[e].isSelected()) return !0;
                    return this._selected;
                },
                setSelected: function (t, e) {
                    if (!e && this._selectChildren)
                        for (
                            var n = this._children, i = 0, r = n.length;
                            r > i;
                            i++
                        )
                            n[i].setSelected(t);
                    (t = !!t) ^ this._selected &&
                        ((this._selected = t),
                        this._project._updateSelection(this),
                        this._changed(129));
                },
                _selected: !1,
                isFullySelected: function () {
                    var t = this._children;
                    if (t && this._selected) {
                        for (var e = 0, n = t.length; n > e; e++)
                            if (!t[e].isFullySelected()) return !1;
                        return !0;
                    }
                    return this._selected;
                },
                setFullySelected: function (t) {
                    var e = this._children;
                    if (e)
                        for (var n = 0, i = e.length; i > n; n++)
                            e[n].setFullySelected(t);
                    this.setSelected(t, !0);
                },
                isClipMask: function () {
                    return this._clipMask;
                },
                setClipMask: function (t) {
                    this._clipMask != (t = !!t) &&
                        ((this._clipMask = t),
                        t &&
                            (this.setFillColor(null),
                            this.setStrokeColor(null)),
                        this._changed(129),
                        this._parent && this._parent._changed(1024));
                },
                _clipMask: !1,
                getData: function () {
                    return this._data || (this._data = {}), this._data;
                },
                setData: function (t) {
                    this._data = t;
                },
                getPosition: function (t) {
                    var e = this._position,
                        n = t ? h : u;
                    if (!e) {
                        var i = this._pivot;
                        e = this._position = i
                            ? this._matrix._transformPoint(i)
                            : this.getBounds().getCenter(!0);
                    }
                    return new n(e.x, e.y, this, "setPosition");
                },
                setPosition: function () {
                    this.translate(
                        h.read(arguments).subtract(this.getPosition(!0))
                    );
                },
                getPivot: function (t) {
                    var e = this._pivot;
                    if (e) {
                        var n = t ? h : u;
                        e = new n(e.x, e.y, this, "setPivot");
                    }
                    return e;
                },
                setPivot: function () {
                    (this._pivot = h.read(arguments)), (this._position = t);
                },
                _pivot: null,
                getRegistration: "#getPivot",
                setRegistration: "#setPivot"
            },
            e.each(
                [
                    "bounds",
                    "strokeBounds",
                    "handleBounds",
                    "roughBounds",
                    "internalBounds",
                    "internalRoughBounds"
                ],
                function (t) {
                    var n = "get" + e.capitalize(t),
                        i = t.match(/^internal(.*)$/),
                        r = i ? "get" + i[1] : null;
                    this[n] = function (e) {
                        var i = this._boundsGetter,
                            s =
                                (!r &&
                                    ("string" == typeof i ? i : i && i[n])) ||
                                n,
                            a = this._getCachedBounds(s, e, this, r);
                        return "bounds" === t
                            ? new _(
                                  a.x,
                                  a.y,
                                  a.width,
                                  a.height,
                                  this,
                                  "setBounds"
                              )
                            : a;
                    };
                },
                {
                    beans: !0,
                    _getBounds: function (t, e, n) {
                        var i = this._children;
                        if (!i || 0 == i.length) return new f();
                        for (
                            var r = 1 / 0,
                                s = -r,
                                a = r,
                                o = s,
                                h = 0,
                                u = i.length;
                            u > h;
                            h++
                        ) {
                            var l = i[h];
                            if (l._visible && !l.isEmpty()) {
                                var c = l._getCachedBounds(t, e, n);
                                (r = Math.min(c.x, r)),
                                    (a = Math.min(c.y, a)),
                                    (s = Math.max(c.x + c.width, s)),
                                    (o = Math.max(c.y + c.height, o));
                            }
                        }
                        return isFinite(r)
                            ? new f(r, a, s - r, o - a)
                            : new f();
                    },
                    setBounds: function () {
                        var t = f.read(arguments),
                            e = this.getBounds(),
                            n = new g(),
                            i = t.getCenter();
                        n.translate(i),
                            (t.width != e.width || t.height != e.height) &&
                                n.scale(
                                    0 != e.width ? t.width / e.width : 1,
                                    0 != e.height ? t.height / e.height : 1
                                ),
                            (i = e.getCenter()),
                            n.translate(-i.x, -i.y),
                            this.transform(n);
                    },
                    _getCachedBounds: function (t, e, n, i) {
                        e = e && e.orNullIfIdentity();
                        var r = i ? null : this._matrix.orNullIfIdentity(),
                            s = (!e || e.equals(r)) && t,
                            a = this._parent || this._parentSymbol;
                        if (a) {
                            var o = n._id,
                                h = (a._boundsCache = a._boundsCache || {
                                    ids: {},
                                    list: []
                                });
                            h.ids[o] || (h.list.push(n), (h.ids[o] = n));
                        }
                        if (s && this._bounds && this._bounds[s])
                            return this._bounds[s].clone();
                        e = e ? (r ? e.chain(r) : e) : r;
                        var u = this._getBounds(i || t, e, n);
                        if (s) {
                            this._bounds || (this._bounds = {});
                            var l = (this._bounds[s] = u.clone());
                            l._internal = !!i;
                        }
                        return u;
                    },
                    statics: {
                        _clearBoundsCache: function (e) {
                            var n = e._boundsCache;
                            if (n) {
                                e._bounds = e._position = e._boundsCache = t;
                                for (
                                    var i = 0, r = n.list, s = r.length;
                                    s > i;
                                    i++
                                ) {
                                    var a = r[i];
                                    a !== e &&
                                        ((a._bounds = a._position = t),
                                        a._boundsCache &&
                                            y._clearBoundsCache(a));
                                }
                            }
                        }
                    }
                }
            ),
            {
                beans: !0,
                _decompose: function () {
                    return (this._decomposed = this._matrix.decompose());
                },
                getRotation: function () {
                    var t = this._decomposed || this._decompose();
                    return t && t.rotation;
                },
                setRotation: function (t) {
                    var e = this.getRotation();
                    if (null != e && null != t) {
                        var n = this._decomposed;
                        this.rotate(t - e),
                            (n.rotation = t),
                            (this._decomposed = n);
                    }
                },
                getScaling: function () {
                    var t = this._decomposed || this._decompose();
                    return t && t.scaling;
                },
                setScaling: function () {
                    var t = this.getScaling();
                    if (null != t) {
                        var e = h.read(arguments, 0, {
                                clone: !0
                            }),
                            n = this._decomposed;
                        this.scale(e.x / t.x, e.y / t.y),
                            (n.scaling = e),
                            (this._decomposed = n);
                    }
                },
                getMatrix: function () {
                    return this._matrix;
                },
                setMatrix: function (t) {
                    this._matrix.initialize(t),
                        this._applyMatrix
                            ? this.transform(null, !0)
                            : this._changed(9);
                },
                getGlobalMatrix: function (t) {
                    var e = this._globalMatrix,
                        n = this._project._updateVersion;
                    if ((e && e._updateVersion !== n && (e = null), !e)) {
                        e = this._globalMatrix = this._matrix.clone();
                        var i = this._parent;
                        i && e.preConcatenate(i.getGlobalMatrix(!0)),
                            (e._updateVersion = n);
                    }
                    return t ? e : e.clone();
                },
                getApplyMatrix: function () {
                    return this._applyMatrix;
                },
                setApplyMatrix: function (t) {
                    (this._applyMatrix = this._canApplyMatrix && !!t) &&
                        this.transform(null, !0);
                },
                getTransformContent: "#getApplyMatrix",
                setTransformContent: "#setApplyMatrix"
            },
            {
                getProject: function () {
                    return this._project;
                },
                _setProject: function (t, e) {
                    if (this._project !== t) {
                        this._project && this._installEvents(!1),
                            (this._project = t);
                        for (
                            var n = this._children, i = 0, r = n && n.length;
                            r > i;
                            i++
                        )
                            n[i]._setProject(t);
                        e = !0;
                    }
                    e && this._installEvents(!0);
                },
                getView: function () {
                    return this._project.getView();
                },
                _installEvents: function le(t) {
                    le.base.call(this, t);
                    for (
                        var e = this._children, n = 0, i = e && e.length;
                        i > n;
                        n++
                    )
                        e[n]._installEvents(t);
                },
                getLayer: function () {
                    for (var t = this; (t = t._parent); )
                        if (t instanceof x) return t;
                    return null;
                },
                getParent: function () {
                    return this._parent;
                },
                setParent: function (t) {
                    return t.addChild(this);
                },
                getChildren: function () {
                    return this._children;
                },
                setChildren: function (t) {
                    this.removeChildren(), this.addChildren(t);
                },
                getFirstChild: function () {
                    return (this._children && this._children[0]) || null;
                },
                getLastChild: function () {
                    return (
                        (this._children &&
                            this._children[this._children.length - 1]) ||
                        null
                    );
                },
                getNextSibling: function () {
                    return (
                        (this._parent &&
                            this._parent._children[this._index + 1]) ||
                        null
                    );
                },
                getPreviousSibling: function () {
                    return (
                        (this._parent &&
                            this._parent._children[this._index - 1]) ||
                        null
                    );
                },
                getIndex: function () {
                    return this._index;
                },
                equals: function (t) {
                    return (
                        t === this ||
                        (t &&
                            this._class === t._class &&
                            this._style.equals(t._style) &&
                            this._matrix.equals(t._matrix) &&
                            this._locked === t._locked &&
                            this._visible === t._visible &&
                            this._blendMode === t._blendMode &&
                            this._opacity === t._opacity &&
                            this._clipMask === t._clipMask &&
                            this._guide === t._guide &&
                            this._equals(t)) ||
                        !1
                    );
                },
                _equals: function (t) {
                    return e.equals(this._children, t._children);
                },
                clone: function (t) {
                    return this._clone(new this.constructor(y.NO_INSERT), t);
                },
                _clone: function (n, i) {
                    if ((n.setStyle(this._style), this._children))
                        for (var r = 0, s = this._children.length; s > r; r++)
                            n.addChild(this._children[r].clone(!1), !0);
                    (i || i === t) && n.insertAbove(this);
                    for (
                        var a = [
                                "_locked",
                                "_visible",
                                "_blendMode",
                                "_opacity",
                                "_clipMask",
                                "_guide",
                                "_applyMatrix"
                            ],
                            r = 0,
                            s = a.length;
                        s > r;
                        r++
                    ) {
                        var o = a[r];
                        this.hasOwnProperty(o) && (n[o] = this[o]);
                    }
                    return (
                        n._matrix.initialize(this._matrix),
                        (n._data = this._data ? e.clone(this._data) : null),
                        n.setSelected(this._selected),
                        this._name && n.setName(this._name, !0),
                        n
                    );
                },
                copyTo: function (t) {
                    return t.addChild(this.clone(!1));
                },
                rasterize: function (t) {
                    var n = this.getStrokeBounds(),
                        i = (t || this.getView().getResolution()) / 72,
                        r = n.getTopLeft().floor(),
                        s = n.getBottomRight().ceil(),
                        a = new c(s.subtract(r)),
                        o = Q.getCanvas(a.multiply(i)),
                        h = o.getContext("2d"),
                        u = new g().scale(i).translate(r.negate());
                    h.save(),
                        u.applyToContext(h),
                        this.draw(
                            h,
                            new e({
                                matrices: [u]
                            })
                        ),
                        h.restore();
                    var l = new C(y.NO_INSERT);
                    return (
                        l.setCanvas(o),
                        l.transform(
                            new g().translate(r.add(a.divide(2))).scale(1 / i)
                        ),
                        l.insertAbove(this),
                        l
                    );
                },
                contains: function () {
                    return !!this._contains(
                        this._matrix._inverseTransform(h.read(arguments))
                    );
                },
                _contains: function (t) {
                    if (this._children) {
                        for (var e = this._children.length - 1; e >= 0; e--)
                            if (this._children[e].contains(t)) return !0;
                        return !1;
                    }
                    return t.isInside(this.getInternalBounds());
                },
                hitTest: function () {
                    return this._hitTest(
                        h.read(arguments),
                        P.getOptions(e.read(arguments))
                    );
                },
                _hitTest: function (n, i) {
                    function r(i, r) {
                        var s = _["get" + r]();
                        return n.subtract(s).divide(u).length <= 1
                            ? new P(i, f, {
                                  name: e.hyphenate(r),
                                  point: s
                              })
                            : t;
                    }
                    if (
                        this._locked ||
                        !this._visible ||
                        (this._guide && !i.guides) ||
                        this.isEmpty()
                    )
                        return null;
                    var s = this._matrix,
                        a = i._totalMatrix,
                        o = this.getView(),
                        h = (i._totalMatrix = a
                            ? a.chain(s)
                            : this.getGlobalMatrix().preConcatenate(o._matrix)),
                        u = (i._tolerancePadding = new c(
                            O._getPenPadding(1, h.inverted())
                        ).multiply(Math.max(i.tolerance, 1e-5)));
                    if (
                        ((n = s._inverseTransform(n)),
                        !this._children &&
                            !this.getInternalRoughBounds()
                                .expand(u.multiply(2))
                                ._containsPoint(n))
                    )
                        return null;
                    var l,
                        d = !(
                            (i.guides && !this._guide) ||
                            (i.selected && !this._selected) ||
                            (i.type && i.type !== e.hyphenate(this._class)) ||
                            (i.class && !(this instanceof i.class))
                        ),
                        f = this;
                    if (d && (i.center || i.bounds) && this._parent) {
                        var _ = this.getInternalBounds();
                        if (
                            (i.center && (l = r("center", "Center")),
                            !l && i.bounds)
                        )
                            for (
                                var g = [
                                        "TopLeft",
                                        "TopRight",
                                        "BottomLeft",
                                        "BottomRight",
                                        "LeftCenter",
                                        "TopCenter",
                                        "RightCenter",
                                        "BottomCenter"
                                    ],
                                    p = 0;
                                8 > p && !l;
                                p++
                            )
                                l = r("bounds", g[p]);
                    }
                    var v = !l && this._children;
                    if (v)
                        for (
                            var m = this._getChildHitTestOptions(i),
                                p = v.length - 1;
                            p >= 0 && !l;
                            p--
                        )
                            l = v[p]._hitTest(n, m);
                    return (
                        !l && d && (l = this._hitTestSelf(n, i)),
                        l && l.point && (l.point = s.transform(l.point)),
                        (i._totalMatrix = a),
                        l
                    );
                },
                _getChildHitTestOptions: function (t) {
                    return t;
                },
                _hitTestSelf: function (e, n) {
                    return n.fill && this.hasFill() && this._contains(e)
                        ? new P("fill", this)
                        : t;
                },
                matches: function (n) {
                    function i(t, n) {
                        for (var r in t)
                            if (t.hasOwnProperty(r)) {
                                var s = t[r],
                                    a = n[r];
                                if (e.isPlainObject(s) && e.isPlainObject(a)) {
                                    if (!i(s, a)) return !1;
                                } else if (!e.equals(s, a)) return !1;
                            }
                        return !0;
                    }
                    for (var r in n)
                        if (n.hasOwnProperty(r)) {
                            var s = this[r],
                                a = n[r];
                            if (
                                (s === t &&
                                    "type" === r &&
                                    (s = e.hyphenate(this._class)),
                                /^(constructor|class)$/.test(r))
                            ) {
                                if (!(this instanceof a)) return !1;
                            } else if (a instanceof RegExp) {
                                if (!a.test(s)) return !1;
                            } else if ("function" == typeof a) {
                                if (!a(s)) return !1;
                            } else if (e.isPlainObject(a)) {
                                if (!i(a, s)) return !1;
                            } else if (!e.equals(s, a)) return !1;
                        }
                    return !0;
                },
                getItems: function (t) {
                    return y._getItems(this._children, t, !0);
                },
                getItem: function (t) {
                    return y._getItems(this._children, t, !1);
                },
                statics: {
                    _getItems: function ce(t, e, n) {
                        for (
                            var i = n && [], r = 0, s = t && t.length;
                            s > r;
                            r++
                        ) {
                            var a = t[r];
                            if (a.matches(e)) {
                                if (!n) return a;
                                i.push(a);
                            }
                            var o = ce(a._children, e, n);
                            if (n) i.push.apply(i, o);
                            else if (o) return o;
                        }
                        return n ? i : null;
                    }
                }
            },
            {
                importJSON: function (t) {
                    var n = e.importJSON(t, this);
                    return n !== this ? this.addChild(n) : n;
                },
                addChild: function (e, n) {
                    return this.insertChild(t, e, n);
                },
                insertChild: function (t, e, n) {
                    var i = this.insertChildren(t, [e], n);
                    return i && i[0];
                },
                addChildren: function (t, e) {
                    return this.insertChildren(this._children.length, t, e);
                },
                insertChildren: function (t, n, i, r) {
                    var s = this._children;
                    if (s && n && n.length > 0) {
                        n = Array.prototype.slice.apply(n);
                        for (var a = n.length - 1; a >= 0; a--) {
                            var o = n[a];
                            !r || o instanceof r
                                ? o._remove(!1, !0)
                                : n.splice(a, 1);
                        }
                        e.splice(s, n, t, 0);
                        for (
                            var h = this._project,
                                u = h && h._changes,
                                a = 0,
                                l = n.length;
                            l > a;
                            a++
                        ) {
                            var o = n[a];
                            (o._parent = this),
                                o._setProject(this._project, !0),
                                o._name && o.setName(o._name),
                                u && this._changed(5);
                        }
                        this._changed(11);
                    } else n = null;
                    return n;
                },
                _insert: function (t, e, n) {
                    if (!e._parent) return null;
                    var i = e._index + (t ? 1 : 0);
                    return (
                        e._parent === this._parent && i > this._index && i--,
                        e._parent.insertChild(i, this, n)
                    );
                },
                insertAbove: function (t, e) {
                    return this._insert(!0, t, e);
                },
                insertBelow: function (t, e) {
                    return this._insert(!1, t, e);
                },
                sendToBack: function () {
                    return this._parent.insertChild(0, this);
                },
                bringToFront: function () {
                    return this._parent.addChild(this);
                },
                appendTop: "#addChild",
                appendBottom: function (t) {
                    return this.insertChild(0, t);
                },
                moveAbove: "#insertAbove",
                moveBelow: "#insertBelow",
                reduce: function () {
                    if (this._children && 1 === this._children.length) {
                        var t = this._children[0].reduce();
                        return (
                            t.insertAbove(this),
                            t.setStyle(this._style),
                            this.remove(),
                            t
                        );
                    }
                    return this;
                },
                _removeNamed: function () {
                    var t = this._parent;
                    if (t) {
                        var e = t._children,
                            n = t._namedChildren,
                            i = this._name,
                            r = n[i],
                            s = r ? r.indexOf(this) : -1;
                        -1 !== s &&
                            (e[i] == this && delete e[i],
                            r.splice(s, 1),
                            r.length ? (e[i] = r[r.length - 1]) : delete n[i]);
                    }
                },
                _remove: function (t, n) {
                    var i = this._parent;
                    if (i) {
                        if (
                            (this._name && this._removeNamed(),
                            null != this._index &&
                                e.splice(i._children, null, this._index, 1),
                            this._installEvents(!1),
                            t)
                        ) {
                            var r = this._project;
                            r && r._changes && this._changed(5);
                        }
                        return n && i._changed(11), (this._parent = null), !0;
                    }
                    return !1;
                },
                remove: function () {
                    return this._remove(!0, !0);
                },
                removeChildren: function (t, n) {
                    if (!this._children) return null;
                    (t = t || 0), (n = e.pick(n, this._children.length));
                    for (
                        var i = e.splice(this._children, null, t, n - t),
                            r = i.length - 1;
                        r >= 0;
                        r--
                    )
                        i[r]._remove(!0, !1);
                    return i.length > 0 && this._changed(11), i;
                },
                clear: "#removeChildren",
                reverseChildren: function () {
                    if (this._children) {
                        this._children.reverse();
                        for (var t = 0, e = this._children.length; e > t; t++)
                            this._children[t]._index = t;
                        this._changed(11);
                    }
                },
                isEmpty: function () {
                    return !this._children || 0 === this._children.length;
                },
                isEditable: function () {
                    for (var t = this; t; ) {
                        if (!t._visible || t._locked) return !1;
                        t = t._parent;
                    }
                    return !0;
                },
                hasFill: function () {
                    return this.getStyle().hasFill();
                },
                hasStroke: function () {
                    return this.getStyle().hasStroke();
                },
                hasShadow: function () {
                    return this.getStyle().hasShadow();
                },
                _getOrder: function (t) {
                    function e(t) {
                        var e = [];
                        do e.unshift(t);
                        while ((t = t._parent));
                        return e;
                    }
                    for (
                        var n = e(this),
                            i = e(t),
                            r = 0,
                            s = Math.min(n.length, i.length);
                        s > r;
                        r++
                    )
                        if (n[r] != i[r])
                            return n[r]._index < i[r]._index ? 1 : -1;
                    return 0;
                },
                hasChildren: function () {
                    return this._children && this._children.length > 0;
                },
                isInserted: function () {
                    return this._parent ? this._parent.isInserted() : !1;
                },
                isAbove: function (t) {
                    return -1 === this._getOrder(t);
                },
                isBelow: function (t) {
                    return 1 === this._getOrder(t);
                },
                isParent: function (t) {
                    return this._parent === t;
                },
                isChild: function (t) {
                    return t && t._parent === this;
                },
                isDescendant: function (t) {
                    for (var e = this; (e = e._parent); ) if (e == t) return !0;
                    return !1;
                },
                isAncestor: function (t) {
                    return t ? t.isDescendant(this) : !1;
                },
                isGroupedWith: function (t) {
                    for (var e = this._parent; e; ) {
                        if (
                            e._parent &&
                            /^(Group|Layer|CompoundPath)$/.test(e._class) &&
                            t.isDescendant(e)
                        )
                            return !0;
                        e = e._parent;
                    }
                    return !1;
                },
                translate: function () {
                    var t = new g();
                    return this.transform(t.translate.apply(t, arguments));
                },
                rotate: function (t) {
                    return this.transform(
                        new g().rotate(
                            t,
                            h.read(arguments, 1, {
                                readNull: !0
                            }) || this.getPosition(!0)
                        )
                    );
                }
            },
            e.each(
                ["scale", "shear", "skew"],
                function (t) {
                    this[t] = function () {
                        var e = h.read(arguments),
                            n = h.read(arguments, 0, {
                                readNull: !0
                            });
                        return this.transform(
                            new g()[t](e, n || this.getPosition(!0))
                        );
                    };
                },
                {}
            ),
            {
                transform: function (t, e) {
                    t && t.isIdentity() && (t = null);
                    var n = this._matrix,
                        i = (e || this._applyMatrix) && (!n.isIdentity() || t);
                    if (!t && !i) return this;
                    if (
                        (t && n.preConcatenate(t),
                        (i = i && this._transformContent(n)))
                    ) {
                        var r = this._pivot,
                            s = this._style,
                            a = s.getFillColor(!0),
                            o = s.getStrokeColor(!0);
                        r && n._transformPoint(r, r, !0),
                            a && a.transform(n),
                            o && o.transform(n),
                            n.reset(!0);
                    }
                    var h = this._bounds,
                        u = this._position;
                    this._changed(9);
                    var l = h && t && t.decompose();
                    if (l && !l.shearing && 0 === l.rotation % 90) {
                        for (var c in h) {
                            var d = h[c];
                            (i || !d._internal) && t._transformBounds(d, d);
                        }
                        var f = this._boundsGetter,
                            d = h[(f && f.getBounds) || f || "getBounds"];
                        d && (this._position = d.getCenter(!0)),
                            (this._bounds = h);
                    } else t && u && (this._position = t._transformPoint(u, u));
                    return this;
                },
                _transformContent: function (t) {
                    var e = this._children;
                    if (e) {
                        for (var n = 0, i = e.length; i > n; n++)
                            e[n].transform(t, !0);
                        return !0;
                    }
                },
                globalToLocal: function () {
                    return this.getGlobalMatrix(!0)._inverseTransform(
                        h.read(arguments)
                    );
                },
                localToGlobal: function () {
                    return this.getGlobalMatrix(!0)._transformPoint(
                        h.read(arguments)
                    );
                },
                fitBounds: function (t, e) {
                    t = f.read(arguments);
                    var n = this.getBounds(),
                        i = n.height / n.width,
                        r = t.height / t.width,
                        s = (e ? i > r : r > i)
                            ? t.width / n.width
                            : t.height / n.height,
                        a = new f(new h(), new c(n.width * s, n.height * s));
                    a.setCenter(t.getCenter()), this.setBounds(a);
                },
                _setStyles: function (t) {
                    var e = this._style,
                        n = e.getFillColor(),
                        i = e.getStrokeColor(),
                        r = e.getShadowColor();
                    if ((n && (t.fillStyle = n.toCanvasStyle(t)), i)) {
                        var s = e.getStrokeWidth();
                        if (s > 0) {
                            (t.strokeStyle = i.toCanvasStyle(t)),
                                (t.lineWidth = s);
                            var a = e.getStrokeJoin(),
                                o = e.getStrokeCap(),
                                h = e.getMiterLimit();
                            if (
                                (a && (t.lineJoin = a),
                                o && (t.lineCap = o),
                                h && (t.miterLimit = h),
                                paper.support.nativeDash)
                            ) {
                                var u = e.getDashArray(),
                                    l = e.getDashOffset();
                                u &&
                                    u.length &&
                                    ("setLineDash" in t
                                        ? (t.setLineDash(u),
                                          (t.lineDashOffset = l))
                                        : ((t.mozDash = u),
                                          (t.mozDashOffset = l)));
                            }
                        }
                    }
                    if (r) {
                        var c = e.getShadowBlur();
                        if (c > 0) {
                            (t.shadowColor = r.toCanvasStyle(t)),
                                (t.shadowBlur = c);
                            var d = this.getShadowOffset();
                            (t.shadowOffsetX = d.x), (t.shadowOffsetY = d.y);
                        }
                    }
                },
                draw: function (t, e, n) {
                    function i(t) {
                        return o ? o.chain(t) : t;
                    }
                    var r = (this._updateVersion =
                        this._project._updateVersion);
                    if (this._visible && 0 !== this._opacity) {
                        var s = e.matrices,
                            a = s[s.length - 1],
                            o = e.viewMatrix,
                            h = this._matrix,
                            u = a.chain(h);
                        if (u.isInvertible()) {
                            s.push(u),
                                e.updateMatrix &&
                                    ((u._updateVersion = r),
                                    (this._globalMatrix = u));
                            var l,
                                c,
                                d,
                                f = this._blendMode,
                                _ = this._opacity,
                                g = "normal" === f,
                                p = te.nativeModes[f],
                                v =
                                    (g && 1 === _) ||
                                    e.dontStart ||
                                    e.clip ||
                                    ((p || (g && 1 > _)) &&
                                        this._canComposite()),
                                m = e.pixelRatio;
                            if (!v) {
                                var y = this.getStrokeBounds(i(a));
                                if (!y.width || !y.height) return;
                                (d = e.offset),
                                    (c = e.offset = y.getTopLeft().floor()),
                                    (l = t),
                                    (t = Q.getContext(
                                        y.getSize().ceil().add(1).multiply(m)
                                    )),
                                    1 !== m && t.scale(m, m);
                            }
                            t.save();
                            var w = n
                                    ? n.chain(h)
                                    : !this.getStrokeScaling(!0) && i(u),
                                x = !v && e.clipItem,
                                b = !w || x;
                            if (
                                (v
                                    ? ((t.globalAlpha = _),
                                      p && (t.globalCompositeOperation = f))
                                    : b && t.translate(-c.x, -c.y),
                                b && (v ? h : i(u)).applyToContext(t),
                                x &&
                                    e.clipItem.draw(
                                        t,
                                        e.extend({
                                            clip: !0
                                        })
                                    ),
                                w)
                            ) {
                                t.setTransform(m, 0, 0, m, 0, 0);
                                var C = e.offset;
                                C && t.translate(-C.x, -C.y);
                            }
                            this._draw(t, e, w),
                                t.restore(),
                                s.pop(),
                                e.clip && !e.dontFinish && t.clip(),
                                v ||
                                    (te.process(
                                        f,
                                        t,
                                        l,
                                        _,
                                        c.subtract(d).multiply(m)
                                    ),
                                    Q.release(t),
                                    (e.offset = d));
                        }
                    }
                },
                _isUpdated: function (t) {
                    var e = this._parent;
                    if (e instanceof T) return e._isUpdated(t);
                    var n = this._updateVersion === t;
                    return (
                        !n &&
                            e &&
                            e._visible &&
                            e._isUpdated(t) &&
                            ((this._updateVersion = t), (n = !0)),
                        n
                    );
                },
                _drawSelection: function (t, e, n, i, r) {
                    if (
                        (this._drawSelected || this._boundsSelected) &&
                        this._isUpdated(r)
                    ) {
                        var s =
                                this.getSelectedColor(!0) ||
                                this.getLayer().getSelectedColor(!0),
                            a = e.chain(this.getGlobalMatrix(!0));
                        if (
                            ((t.strokeStyle = t.fillStyle =
                                s ? s.toCanvasStyle(t) : "#009dec"),
                            this._drawSelected && this._drawSelected(t, a, i),
                            this._boundsSelected)
                        ) {
                            var o = n / 2;
                            (coords = a._transformCorners(
                                this.getInternalBounds()
                            )),
                                t.beginPath();
                            for (var h = 0; 8 > h; h++)
                                t[0 === h ? "moveTo" : "lineTo"](
                                    coords[h],
                                    coords[++h]
                                );
                            t.closePath(), t.stroke();
                            for (var h = 0; 8 > h; h++)
                                t.fillRect(
                                    coords[h] - o,
                                    coords[++h] - o,
                                    n,
                                    n
                                );
                        }
                    }
                },
                _canComposite: function () {
                    return !1;
                }
            },
            e.each(
                ["down", "drag", "up", "move"],
                function (t) {
                    this["removeOn" + e.capitalize(t)] = function () {
                        var e = {};
                        return (e[t] = !0), this.removeOn(e);
                    };
                },
                {
                    removeOn: function (t) {
                        for (var e in t)
                            if (t[e]) {
                                var n = "mouse" + e,
                                    i = this._project,
                                    r = (i._removeSets = i._removeSets || {});
                                (r[n] = r[n] || {}), (r[n][this._id] = this);
                            }
                        return this;
                    }
                }
            )
        ),
        w = y.extend({
            _class: "Group",
            _selectChildren: !0,
            _serializeFields: {
                children: []
            },
            initialize: function (t) {
                (this._children = []),
                    (this._namedChildren = {}),
                    this._initialize(t) ||
                        this.addChildren(Array.isArray(t) ? t : arguments);
            },
            _changed: function de(e) {
                de.base.call(this, e), 1026 & e && (this._clipItem = t);
            },
            _getClipItem: function () {
                var e = this._clipItem;
                if (e === t) {
                    e = null;
                    for (var n = 0, i = this._children.length; i > n; n++) {
                        var r = this._children[n];
                        if (r._clipMask) {
                            e = r;
                            break;
                        }
                    }
                    this._clipItem = e;
                }
                return e;
            },
            isClipped: function () {
                return !!this._getClipItem();
            },
            setClipped: function (t) {
                var e = this.getFirstChild();
                e && e.setClipMask(t);
            },
            _draw: function (t, e) {
                var n = e.clip,
                    i = !n && this._getClipItem(),
                    r = !0;
                if (
                    ((e = e.extend({
                        clipItem: i,
                        clip: !1
                    })),
                    n
                        ? this._currentPath
                            ? ((t.currentPath = this._currentPath), (r = !1))
                            : (t.beginPath(), (e.dontStart = e.dontFinish = !0))
                        : i &&
                          i.draw(
                              t,
                              e.extend({
                                  clip: !0
                              })
                          ),
                    r)
                )
                    for (var s = 0, a = this._children.length; a > s; s++) {
                        var o = this._children[s];
                        o !== i && o.draw(t, e);
                    }
                n && (this._currentPath = t.currentPath);
            }
        }),
        x = w.extend({
            _class: "Layer",
            initialize: function (n) {
                var i = e.isPlainObject(n)
                        ? new e(n)
                        : {
                              children: Array.isArray(n) ? n : arguments
                          },
                    r = i.insert;
                (i.insert = !1),
                    w.call(this, i),
                    (r || r === t) &&
                        (this._project.addChild(this), this.activate());
            },
            _remove: function fe(t) {
                return this._parent
                    ? fe.base.call(this, t)
                    : null != this._index
                    ? (this._project.activeLayer === this &&
                          (this._project.activeLayer =
                              this.getNextSibling() ||
                              this.getPreviousSibling()),
                      e.splice(this._project.layers, null, this._index, 1),
                      this._installEvents(!1),
                      (this._project._needsUpdate = !0),
                      !0)
                    : !1;
            },
            getNextSibling: function _e() {
                return this._parent
                    ? _e.base.call(this)
                    : this._project.layers[this._index + 1] || null;
            },
            getPreviousSibling: function ge() {
                return this._parent
                    ? ge.base.call(this)
                    : this._project.layers[this._index - 1] || null;
            },
            isInserted: function pe() {
                return this._parent ? pe.base.call(this) : null != this._index;
            },
            activate: function () {
                this._project.activeLayer = this;
            },
            _insert: function ve(t, n, i) {
                return n instanceof x && !n._parent
                    ? (this._remove(!0, !0),
                      e.splice(
                          n._project.layers,
                          [this],
                          n._index + (t ? 1 : 0),
                          0
                      ),
                      this._setProject(n._project, !0),
                      this)
                    : ve.base.call(this, t, n, i);
            }
        }),
        b = y.extend(
            {
                _class: "Shape",
                _applyMatrix: !1,
                _canApplyMatrix: !1,
                _boundsSelected: !0,
                _serializeFields: {
                    type: null,
                    size: null,
                    radius: null
                },
                initialize: function (t) {
                    this._initialize(t);
                },
                _equals: function (t) {
                    return (
                        this._type === t._type &&
                        this._size.equals(t._size) &&
                        e.equals(this._radius, t._radius)
                    );
                },
                clone: function (t) {
                    var e = new b(y.NO_INSERT);
                    return (
                        e.setType(this._type),
                        e.setSize(this._size),
                        e.setRadius(this._radius),
                        this._clone(e, t)
                    );
                },
                getType: function () {
                    return this._type;
                },
                setType: function (t) {
                    this._type = t;
                },
                getShape: "#getType",
                setShape: "#setType",
                getSize: function () {
                    var t = this._size;
                    return new d(t.width, t.height, this, "setSize");
                },
                setSize: function () {
                    var t = c.read(arguments);
                    if (this._size) {
                        if (!this._size.equals(t)) {
                            var e = this._type,
                                n = t.width,
                                i = t.height;
                            if ("rectangle" === e) {
                                var r = c.min(this._radius, t.divide(2));
                                this._radius.set(r.width, r.height);
                            } else
                                "circle" === e
                                    ? ((n = i = (n + i) / 2),
                                      (this._radius = n / 2))
                                    : "ellipse" === e &&
                                      this._radius.set(n / 2, i / 2);
                            this._size.set(n, i), this._changed(9);
                        }
                    } else this._size = t.clone();
                },
                getRadius: function () {
                    var t = this._radius;
                    return "circle" === this._type
                        ? t
                        : new d(t.width, t.height, this, "setRadius");
                },
                setRadius: function (t) {
                    var e = this._type;
                    if ("circle" === e) {
                        if (t === this._radius) return;
                        var n = 2 * t;
                        (this._radius = t), this._size.set(n, n);
                    } else if (((t = c.read(arguments)), this._radius)) {
                        if (this._radius.equals(t)) return;
                        if (
                            (this._radius.set(t.width, t.height),
                            "rectangle" === e)
                        ) {
                            var n = c.max(this._size, t.multiply(2));
                            this._size.set(n.width, n.height);
                        } else
                            "ellipse" === e &&
                                this._size.set(2 * t.width, 2 * t.height);
                    } else this._radius = t.clone();
                    this._changed(9);
                },
                isEmpty: function () {
                    return !1;
                },
                toPath: function (n) {
                    var i = new O[e.capitalize(this._type)]({
                        center: new h(),
                        size: this._size,
                        radius: this._radius,
                        insert: !1
                    });
                    return (
                        i.setStyle(this._style),
                        i.transform(this._matrix),
                        (n || n === t) && i.insertAbove(this),
                        i
                    );
                },
                _draw: function (t, e, n) {
                    var i = this._style,
                        r = i.hasFill(),
                        s = i.hasStroke(),
                        a = e.dontFinish || e.clip,
                        o = !n;
                    if (r || s || a) {
                        var h = this._type,
                            u = this._radius,
                            l = "circle" === h;
                        if ((e.dontStart || t.beginPath(), o && l))
                            t.arc(0, 0, u, 0, 2 * Math.PI, !0);
                        else {
                            var c = l ? u : u.width,
                                d = l ? u : u.height,
                                f = this._size,
                                _ = f.width,
                                g = f.height;
                            if (o && "rect" === h && 0 === c && 0 === d)
                                t.rect(-_ / 2, -g / 2, _, g);
                            else {
                                var p = _ / 2,
                                    v = g / 2,
                                    m = 0.44771525016920644,
                                    y = c * m,
                                    w = d * m,
                                    x = [
                                        -p,
                                        -v + d,
                                        -p,
                                        -v + w,
                                        -p + y,
                                        -v,
                                        -p + c,
                                        -v,
                                        p - c,
                                        -v,
                                        p - y,
                                        -v,
                                        p,
                                        -v + w,
                                        p,
                                        -v + d,
                                        p,
                                        v - d,
                                        p,
                                        v - w,
                                        p - y,
                                        v,
                                        p - c,
                                        v,
                                        -p + c,
                                        v,
                                        -p + y,
                                        v,
                                        -p,
                                        v - w,
                                        -p,
                                        v - d
                                    ];
                                n && n.transform(x, x, 32),
                                    t.moveTo(x[0], x[1]),
                                    t.bezierCurveTo(
                                        x[2],
                                        x[3],
                                        x[4],
                                        x[5],
                                        x[6],
                                        x[7]
                                    ),
                                    p !== c && t.lineTo(x[8], x[9]),
                                    t.bezierCurveTo(
                                        x[10],
                                        x[11],
                                        x[12],
                                        x[13],
                                        x[14],
                                        x[15]
                                    ),
                                    v !== d && t.lineTo(x[16], x[17]),
                                    t.bezierCurveTo(
                                        x[18],
                                        x[19],
                                        x[20],
                                        x[21],
                                        x[22],
                                        x[23]
                                    ),
                                    p !== c && t.lineTo(x[24], x[25]),
                                    t.bezierCurveTo(
                                        x[26],
                                        x[27],
                                        x[28],
                                        x[29],
                                        x[30],
                                        x[31]
                                    );
                            }
                        }
                        t.closePath();
                    }
                    a ||
                        (!r && !s) ||
                        (this._setStyles(t),
                        r &&
                            (t.fill(i.getWindingRule()),
                            (t.shadowColor = "rgba(0,0,0,0)")),
                        s && t.stroke());
                },
                _canComposite: function () {
                    return !(this.hasFill() && this.hasStroke());
                },
                _getBounds: function (t, e) {
                    var n = new f(this._size).setCenter(0, 0);
                    return (
                        "getBounds" !== t &&
                            this.hasStroke() &&
                            (n = n.expand(this.getStrokeWidth())),
                        e ? e._transformBounds(n) : n
                    );
                }
            },
            new (function () {
                function t(t, e, n) {
                    var i = t._radius;
                    if (!i.isZero())
                        for (var r = t._size.divide(2), s = 0; 4 > s; s++) {
                            var a = new h(1 & s ? 1 : -1, s > 1 ? 1 : -1),
                                o = a.multiply(r),
                                u = o.subtract(a.multiply(i)),
                                l = new f(o, u);
                            if ((n ? l.expand(n) : l).contains(e)) return u;
                        }
                }
                function e(t, e) {
                    var n = t.getAngleInRadians(),
                        i = 2 * e.width,
                        r = 2 * e.height,
                        s = i * Math.sin(n),
                        a = r * Math.cos(n);
                    return (i * r) / (2 * Math.sqrt(s * s + a * a));
                }
                return {
                    _contains: function n(e) {
                        if ("rectangle" === this._type) {
                            var i = t(this, e);
                            return i
                                ? e
                                      .subtract(i)
                                      .divide(this._radius)
                                      .getLength() <= 1
                                : n.base.call(this, e);
                        }
                        return e.divide(this.size).getLength() <= 0.5;
                    },
                    _hitTestSelf: function i(n, r) {
                        var s = !1;
                        if (this.hasStroke()) {
                            var a = this._type,
                                o = this._radius,
                                h = this.getStrokeWidth() + 2 * r.tolerance;
                            if ("rectangle" === a) {
                                var u = t(this, n, h);
                                if (u) {
                                    var l = n.subtract(u);
                                    s =
                                        2 * Math.abs(l.getLength() - e(l, o)) <=
                                        h;
                                } else {
                                    var c = new f(this._size).setCenter(0, 0),
                                        d = c.expand(h),
                                        _ = c.expand(-h);
                                    s =
                                        d._containsPoint(n) &&
                                        !_._containsPoint(n);
                                }
                            } else
                                "ellipse" === a && (o = e(n, o)),
                                    (s = 2 * Math.abs(n.getLength() - o) <= h);
                        }
                        return s
                            ? new P("stroke", this)
                            : i.base.apply(this, arguments);
                    }
                };
            })(),
            {
                statics: new (function () {
                    function t(t, n, i, r, s) {
                        var a = new b(e.getNamed(s));
                        return (
                            (a._type = t),
                            (a._size = i),
                            (a._radius = r),
                            a.translate(n)
                        );
                    }
                    return {
                        Circle: function () {
                            var n = h.readNamed(arguments, "center"),
                                i = e.readNamed(arguments, "radius");
                            return t("circle", n, new c(2 * i), i, arguments);
                        },
                        Rectangle: function () {
                            var e = f.readNamed(arguments, "rectangle"),
                                n = c.min(
                                    c.readNamed(arguments, "radius"),
                                    e.getSize(!0).divide(2)
                                );
                            return t(
                                "rectangle",
                                e.getCenter(!0),
                                e.getSize(!0),
                                n,
                                arguments
                            );
                        },
                        Ellipse: function () {
                            var e = b._readEllipse(arguments),
                                n = e.radius;
                            return t(
                                "ellipse",
                                e.center,
                                n.multiply(2),
                                n,
                                arguments
                            );
                        },
                        _readEllipse: function (t) {
                            var n, i;
                            if (e.hasNamed(t, "radius"))
                                (n = h.readNamed(t, "center")),
                                    (i = c.readNamed(t, "radius"));
                            else {
                                var r = f.readNamed(t, "rectangle");
                                (n = r.getCenter(!0)),
                                    (i = r.getSize(!0).divide(2));
                            }
                            return {
                                center: n,
                                radius: i
                            };
                        }
                    };
                })()
            }
        ),
        C = y.extend(
            {
                _class: "Raster",
                _applyMatrix: !1,
                _canApplyMatrix: !1,
                _boundsGetter: "getBounds",
                _boundsSelected: !0,
                _serializeFields: {
                    source: null
                },
                initialize: function (e, n) {
                    this._initialize(e, n !== t && h.read(arguments, 1)) ||
                        ("string" == typeof e
                            ? this.setSource(e)
                            : this.setImage(e)),
                        this._size || (this._size = new c());
                },
                _equals: function (t) {
                    return this.getSource() === t.getSource();
                },
                clone: function (t) {
                    var e = new C(y.NO_INSERT),
                        n = this._image,
                        i = this._canvas;
                    if (n) e.setImage(n);
                    else if (i) {
                        var r = Q.getCanvas(this._size);
                        r.getContext("2d").drawImage(i, 0, 0), e.setCanvas(r);
                    }
                    return this._clone(e, t);
                },
                getSize: function () {
                    var t = this._size;
                    return new d(t.width, t.height, this, "setSize");
                },
                setSize: function () {
                    var t = c.read(arguments);
                    if (!this._size.equals(t)) {
                        var e = this.getElement();
                        this.setCanvas(Q.getCanvas(t)),
                            e &&
                                this.getContext(!0).drawImage(
                                    e,
                                    0,
                                    0,
                                    t.width,
                                    t.height
                                );
                    }
                },
                getWidth: function () {
                    return this._size.width;
                },
                getHeight: function () {
                    return this._size.height;
                },
                isEmpty: function () {
                    return 0 === this._size.width && 0 === this._size.height;
                },
                getResolution: function () {
                    var t = this._matrix,
                        e = new h(0, 0).transform(t),
                        n = new h(1, 0).transform(t).subtract(e),
                        i = new h(0, 1).transform(t).subtract(e);
                    return new c(72 / n.getLength(), 72 / i.getLength());
                },
                getPpi: "#getResolution",
                getImage: function () {
                    return this._image;
                },
                setImage: function (t) {
                    this._canvas && Q.release(this._canvas),
                        t && t.getContext
                            ? ((this._image = null), (this._canvas = t))
                            : ((this._image = t), (this._canvas = null)),
                        (this._size = new c(
                            t ? t.naturalWidth || t.width : 0,
                            t ? t.naturalHeight || t.height : 0
                        )),
                        (this._context = null),
                        this._changed(521);
                },
                getCanvas: function () {
                    if (!this._canvas) {
                        var t = Q.getContext(this._size);
                        try {
                            this._image && t.drawImage(this._image, 0, 0),
                                (this._canvas = t.canvas);
                        } catch (e) {
                            Q.release(t);
                        }
                    }
                    return this._canvas;
                },
                setCanvas: "#setImage",
                getContext: function (t) {
                    return (
                        this._context ||
                            (this._context = this.getCanvas().getContext("2d")),
                        t && ((this._image = null), this._changed(513)),
                        this._context
                    );
                },
                setContext: function (t) {
                    this._context = t;
                },
                getSource: function () {
                    return (this._image && this._image.src) || this.toDataURL();
                },
                setSource: function (t) {
                    function e() {
                        var t = i.getView();
                        t &&
                            ((paper = t._scope),
                            i.setImage(n),
                            i.fire("load"),
                            t.update());
                    }
                    var n,
                        i = this;
                    (n = document.getElementById(t) || new Image()),
                        n.naturalWidth && n.naturalHeight
                            ? setTimeout(e, 0)
                            : (q.add(n, {
                                  load: e
                              }),
                              n.src || (n.src = t)),
                        this.setImage(n);
                },
                getElement: function () {
                    return this._canvas || this._image;
                }
            },
            {
                beans: !1,
                getSubCanvas: function () {
                    var t = f.read(arguments),
                        e = Q.getContext(t.getSize());
                    return (
                        e.drawImage(
                            this.getCanvas(),
                            t.x,
                            t.y,
                            t.width,
                            t.height,
                            0,
                            0,
                            t.width,
                            t.height
                        ),
                        e.canvas
                    );
                },
                getSubRaster: function () {
                    var t = f.read(arguments),
                        e = new C(y.NO_INSERT);
                    return (
                        e.setCanvas(this.getSubCanvas(t)),
                        e.translate(
                            t.getCenter().subtract(this.getSize().divide(2))
                        ),
                        e._matrix.preConcatenate(this._matrix),
                        e.insertAbove(this),
                        e
                    );
                },
                toDataURL: function () {
                    var t = this._image && this._image.src;
                    if (/^data:/.test(t)) return t;
                    var e = this.getCanvas();
                    return e ? e.toDataURL() : null;
                },
                drawImage: function (t) {
                    var e = h.read(arguments, 1);
                    this.getContext(!0).drawImage(t, e.x, e.y);
                },
                getAverageColor: function (t) {
                    var n, i;
                    t
                        ? t instanceof I
                            ? ((i = t), (n = t.getBounds()))
                            : t.width
                            ? (n = new f(t))
                            : t.x && (n = new f(t.x - 0.5, t.y - 0.5, 1, 1))
                        : (n = this.getBounds());
                    var r = 32,
                        s = Math.min(n.width, r),
                        a = Math.min(n.height, r),
                        o = C._sampleContext;
                    o
                        ? o.clearRect(0, 0, r + 1, r + 1)
                        : (o = C._sampleContext = Q.getContext(new c(r))),
                        o.save();
                    var h = new g()
                        .scale(s / n.width, a / n.height)
                        .translate(-n.x, -n.y);
                    h.applyToContext(o),
                        i &&
                            i.draw(
                                o,
                                new e({
                                    clip: !0,
                                    matrices: [h]
                                })
                            ),
                        this._matrix.applyToContext(o),
                        o.drawImage(
                            this.getElement(),
                            -this._size.width / 2,
                            -this._size.height / 2
                        ),
                        o.restore();
                    for (
                        var u = o.getImageData(
                                0.5,
                                0.5,
                                Math.ceil(s),
                                Math.ceil(a)
                            ).data,
                            l = [0, 0, 0],
                            d = 0,
                            _ = 0,
                            p = u.length;
                        p > _;
                        _ += 4
                    ) {
                        var v = u[_ + 3];
                        (d += v),
                            (v /= 255),
                            (l[0] += u[_] * v),
                            (l[1] += u[_ + 1] * v),
                            (l[2] += u[_ + 2] * v);
                    }
                    for (var _ = 0; 3 > _; _++) l[_] /= d;
                    return d ? D.read(l) : null;
                },
                getPixel: function () {
                    var t = h.read(arguments),
                        e = this.getContext().getImageData(t.x, t.y, 1, 1).data;
                    return new D(
                        "rgb",
                        [e[0] / 255, e[1] / 255, e[2] / 255],
                        e[3] / 255
                    );
                },
                setPixel: function () {
                    var t = h.read(arguments),
                        e = D.read(arguments),
                        n = e._convert("rgb"),
                        i = e._alpha,
                        r = this.getContext(!0),
                        s = r.createImageData(1, 1),
                        a = s.data;
                    (a[0] = 255 * n[0]),
                        (a[1] = 255 * n[1]),
                        (a[2] = 255 * n[2]),
                        (a[3] = null != i ? 255 * i : 255),
                        r.putImageData(s, t.x, t.y);
                },
                createImageData: function () {
                    var t = c.read(arguments);
                    return this.getContext().createImageData(t.width, t.height);
                },
                getImageData: function () {
                    var t = f.read(arguments);
                    return (
                        t.isEmpty() && (t = new f(this._size)),
                        this.getContext().getImageData(
                            t.x,
                            t.y,
                            t.width,
                            t.height
                        )
                    );
                },
                setImageData: function (t) {
                    var e = h.read(arguments, 1);
                    this.getContext(!0).putImageData(t, e.x, e.y);
                },
                _getBounds: function (t, e) {
                    var n = new f(this._size).setCenter(0, 0);
                    return e ? e._transformBounds(n) : n;
                },
                _hitTestSelf: function (t) {
                    if (this._contains(t)) {
                        var e = this;
                        return new P("pixel", e, {
                            offset: t.add(e._size.divide(2)).round(),
                            color: {
                                get: function () {
                                    return e.getPixel(this.offset);
                                }
                            }
                        });
                    }
                },
                _draw: function (t) {
                    var e = this.getElement();
                    e &&
                        ((t.globalAlpha = this._opacity),
                        t.drawImage(
                            e,
                            -this._size.width / 2,
                            -this._size.height / 2
                        ));
                },
                _canComposite: function () {
                    return !0;
                }
            }
        ),
        S = y.extend({
            _class: "PlacedSymbol",
            _applyMatrix: !1,
            _canApplyMatrix: !1,
            _boundsGetter: {
                getBounds: "getStrokeBounds"
            },
            _boundsSelected: !0,
            _serializeFields: {
                symbol: null
            },
            initialize: function (e, n) {
                this._initialize(e, n !== t && h.read(arguments, 1)) ||
                    this.setSymbol(e instanceof m ? e : new m(e));
            },
            _equals: function (t) {
                return this._symbol === t._symbol;
            },
            getSymbol: function () {
                return this._symbol;
            },
            setSymbol: function (t) {
                (this._symbol = t), this._changed(9);
            },
            clone: function (t) {
                var e = new S(y.NO_INSERT);
                return e.setSymbol(this._symbol), this._clone(e, t);
            },
            isEmpty: function () {
                return this._symbol._definition.isEmpty();
            },
            _getBounds: function (t, e, n) {
                return this.symbol._definition._getCachedBounds(t, e, n);
            },
            _hitTestSelf: function (t, e) {
                var n = this._symbol._definition._hitTest(t, e);
                return n && (n.item = this), n;
            },
            _draw: function (t, e) {
                this.symbol._definition.draw(t, e);
            }
        }),
        P = e.extend({
            _class: "HitResult",
            initialize: function (t, e, n) {
                (this.type = t),
                    (this.item = e),
                    n && ((n.enumerable = !0), this.inject(n));
            },
            statics: {
                getOptions: function (t) {
                    return new e(
                        {
                            type: null,
                            tolerance: paper.settings.hitTolerance,
                            fill: !t,
                            stroke: !t,
                            segments: !t,
                            handles: !1,
                            ends: !1,
                            center: !1,
                            bounds: !1,
                            guides: !1,
                            selected: !1
                        },
                        t
                    );
                }
            }
        }),
        k = e.extend({
            _class: "Segment",
            beans: !0,
            initialize: function (e, n, i, r, s, a) {
                var o,
                    h,
                    u,
                    l = arguments.length;
                0 === l ||
                    (1 === l
                        ? e.point
                            ? ((o = e.point),
                              (h = e.handleIn),
                              (u = e.handleOut))
                            : (o = e)
                        : 2 === l && "number" == typeof e
                        ? (o = arguments)
                        : 3 >= l
                        ? ((o = e), (h = n), (u = i))
                        : ((o = e !== t ? [e, n] : null),
                          (h = i !== t ? [i, r] : null),
                          (u = s !== t ? [s, a] : null))),
                    new M(o, this, "_point"),
                    new M(h, this, "_handleIn"),
                    new M(u, this, "_handleOut");
            },
            _serialize: function (t) {
                return e.serialize(
                    this.isLinear()
                        ? this._point
                        : [this._point, this._handleIn, this._handleOut],
                    t,
                    !0
                );
            },
            _changed: function (t) {
                var e = this._path;
                if (e) {
                    var n,
                        i,
                        r = e._curves,
                        s = this._index;
                    r &&
                        ((t && t !== this._point && t !== this._handleIn) ||
                            !(n = r[s - 1] || (e._closed && r[r.length - 1])) ||
                            n._changed(),
                        (t && t !== this._point && t !== this._handleOut) ||
                            !(i = r[s]) ||
                            i._changed()),
                        e._changed(25);
                }
            },
            getPoint: function () {
                return this._point;
            },
            setPoint: function () {
                var t = h.read(arguments);
                this._point.set(t.x, t.y);
            },
            getHandleIn: function () {
                return this._handleIn;
            },
            setHandleIn: function () {
                var t = h.read(arguments);
                this._handleIn.set(t.x, t.y);
            },
            getHandleOut: function () {
                return this._handleOut;
            },
            setHandleOut: function () {
                var t = h.read(arguments);
                this._handleOut.set(t.x, t.y);
            },
            isLinear: function () {
                return this._handleIn.isZero() && this._handleOut.isZero();
            },
            setLinear: function (t) {
                t && (this._handleIn.set(0, 0), this._handleOut.set(0, 0));
            },
            isColinear: function (t) {
                var e = this.getNext(),
                    n = t.getNext();
                return (
                    this._handleOut.isZero() &&
                    e._handleIn.isZero() &&
                    t._handleOut.isZero() &&
                    n._handleIn.isZero() &&
                    e._point
                        .subtract(this._point)
                        .isColinear(n._point.subtract(t._point))
                );
            },
            isOrthogonal: function () {
                var t = this.getPrevious(),
                    e = this.getNext();
                return (
                    t._handleOut.isZero() &&
                    this._handleIn.isZero() &&
                    this._handleOut.isZero() &&
                    e._handleIn.isZero() &&
                    this._point
                        .subtract(t._point)
                        .isOrthogonal(e._point.subtract(this._point))
                );
            },
            isArc: function () {
                var t = this.getNext(),
                    e = this._handleOut,
                    n = t._handleIn,
                    i = 0.5522847498307936;
                if (e.isOrthogonal(n)) {
                    var r = this._point,
                        s = t._point,
                        a = new p(r, e, !0).intersect(new p(s, n, !0), !0);
                    return (
                        a &&
                        o.isZero(
                            e.getLength() / a.subtract(r).getLength() - i
                        ) &&
                        o.isZero(n.getLength() / a.subtract(s).getLength() - i)
                    );
                }
                return !1;
            },
            _selectionState: 0,
            isSelected: function (t) {
                var e = this._selectionState;
                return t
                    ? t === this._point
                        ? !!(4 & e)
                        : t === this._handleIn
                        ? !!(1 & e)
                        : t === this._handleOut
                        ? !!(2 & e)
                        : !1
                    : !!(7 & e);
            },
            setSelected: function (t, e) {
                var n = this._path,
                    t = !!t,
                    i = this._selectionState,
                    r = i,
                    s = e
                        ? e === this._point
                            ? 4
                            : e === this._handleIn
                            ? 1
                            : e === this._handleOut
                            ? 2
                            : 0
                        : 7;
                t ? (i |= s) : (i &= ~s),
                    (this._selectionState = i),
                    n &&
                        i !== r &&
                        (n._updateSelection(this, r, i), n._changed(129));
            },
            getIndex: function () {
                return this._index !== t ? this._index : null;
            },
            getPath: function () {
                return this._path || null;
            },
            getCurve: function () {
                var t = this._path,
                    e = this._index;
                return t
                    ? (e > 0 &&
                          !t._closed &&
                          e === t._segments.length - 1 &&
                          e--,
                      t.getCurves()[e] || null)
                    : null;
            },
            getLocation: function () {
                var t = this.getCurve();
                return t ? new A(t, this === t._segment1 ? 0 : 1) : null;
            },
            getNext: function () {
                var t = this._path && this._path._segments;
                return (
                    (t &&
                        (t[this._index + 1] || (this._path._closed && t[0]))) ||
                    null
                );
            },
            getPrevious: function () {
                var t = this._path && this._path._segments;
                return (
                    (t &&
                        (t[this._index - 1] ||
                            (this._path._closed && t[t.length - 1]))) ||
                    null
                );
            },
            reverse: function () {
                return new k(this._point, this._handleOut, this._handleIn);
            },
            remove: function () {
                return this._path
                    ? !!this._path.removeSegment(this._index)
                    : !1;
            },
            clone: function () {
                return new k(this._point, this._handleIn, this._handleOut);
            },
            equals: function (t) {
                return (
                    t === this ||
                    (t &&
                        this._class === t._class &&
                        this._point.equals(t._point) &&
                        this._handleIn.equals(t._handleIn) &&
                        this._handleOut.equals(t._handleOut)) ||
                    !1
                );
            },
            toString: function () {
                var t = ["point: " + this._point];
                return (
                    this._handleIn.isZero() ||
                        t.push("handleIn: " + this._handleIn),
                    this._handleOut.isZero() ||
                        t.push("handleOut: " + this._handleOut),
                    "{ " + t.join(", ") + " }"
                );
            },
            transform: function (t) {
                this._transformCoordinates(t, Array(6), !0), this._changed();
            },
            _transformCoordinates: function (t, e, n) {
                var i = this._point,
                    r = n && this._handleIn.isZero() ? null : this._handleIn,
                    s = n && this._handleOut.isZero() ? null : this._handleOut,
                    a = i._x,
                    o = i._y,
                    h = 2;
                return (
                    (e[0] = a),
                    (e[1] = o),
                    r && ((e[h++] = r._x + a), (e[h++] = r._y + o)),
                    s && ((e[h++] = s._x + a), (e[h++] = s._y + o)),
                    t &&
                        (t._transformCoordinates(e, e, h / 2),
                        (a = e[0]),
                        (o = e[1]),
                        n
                            ? ((i._x = a),
                              (i._y = o),
                              (h = 2),
                              r && ((r._x = e[h++] - a), (r._y = e[h++] - o)),
                              s && ((s._x = e[h++] - a), (s._y = e[h++] - o)))
                            : (r || ((e[h++] = a), (e[h++] = o)),
                              s || ((e[h++] = a), (e[h++] = o)))),
                    e
                );
            }
        }),
        M = h.extend({
            initialize: function (e, n, i) {
                var r, s, a;
                if (e)
                    if ((r = e[0]) !== t) s = e[1];
                    else {
                        var o = e;
                        (r = o.x) === t && ((o = h.read(arguments)), (r = o.x)),
                            (s = o.y),
                            (a = o.selected);
                    }
                else r = s = 0;
                (this._x = r),
                    (this._y = s),
                    (this._owner = n),
                    (n[i] = this),
                    a && this.setSelected(!0);
            },
            set: function (t, e) {
                return (
                    (this._x = t),
                    (this._y = e),
                    this._owner._changed(this),
                    this
                );
            },
            _serialize: function (t) {
                var e = t.formatter,
                    n = e.number(this._x),
                    i = e.number(this._y);
                return this.isSelected()
                    ? {
                          x: n,
                          y: i,
                          selected: !0
                      }
                    : [n, i];
            },
            getX: function () {
                return this._x;
            },
            setX: function (t) {
                (this._x = t), this._owner._changed(this);
            },
            getY: function () {
                return this._y;
            },
            setY: function (t) {
                (this._y = t), this._owner._changed(this);
            },
            isZero: function () {
                return o.isZero(this._x) && o.isZero(this._y);
            },
            setSelected: function (t) {
                this._owner.setSelected(t, this);
            },
            isSelected: function () {
                return this._owner.isSelected(this);
            }
        }),
        z = e.extend(
            {
                _class: "Curve",
                initialize: function (t, e, n, i, r, s, a, o) {
                    var h = arguments.length;
                    if (3 === h)
                        (this._path = t),
                            (this._segment1 = e),
                            (this._segment2 = n);
                    else if (0 === h)
                        (this._segment1 = new k()), (this._segment2 = new k());
                    else if (1 === h)
                        (this._segment1 = new k(t.segment1)),
                            (this._segment2 = new k(t.segment2));
                    else if (2 === h)
                        (this._segment1 = new k(t)),
                            (this._segment2 = new k(e));
                    else {
                        var u, l, c, d;
                        4 === h
                            ? ((u = t), (l = e), (c = n), (d = i))
                            : 8 === h &&
                              ((u = [t, e]),
                              (d = [a, o]),
                              (l = [n - t, i - e]),
                              (c = [r - a, s - o])),
                            (this._segment1 = new k(u, null, l)),
                            (this._segment2 = new k(d, c, null));
                    }
                },
                _changed: function () {
                    this._length = this._bounds = t;
                },
                getPoint1: function () {
                    return this._segment1._point;
                },
                setPoint1: function () {
                    var t = h.read(arguments);
                    this._segment1._point.set(t.x, t.y);
                },
                getPoint2: function () {
                    return this._segment2._point;
                },
                setPoint2: function () {
                    var t = h.read(arguments);
                    this._segment2._point.set(t.x, t.y);
                },
                getHandle1: function () {
                    return this._segment1._handleOut;
                },
                setHandle1: function () {
                    var t = h.read(arguments);
                    this._segment1._handleOut.set(t.x, t.y);
                },
                getHandle2: function () {
                    return this._segment2._handleIn;
                },
                setHandle2: function () {
                    var t = h.read(arguments);
                    this._segment2._handleIn.set(t.x, t.y);
                },
                getSegment1: function () {
                    return this._segment1;
                },
                getSegment2: function () {
                    return this._segment2;
                },
                getPath: function () {
                    return this._path;
                },
                getIndex: function () {
                    return this._segment1._index;
                },
                getNext: function () {
                    var t = this._path && this._path._curves;
                    return (
                        (t &&
                            (t[this._segment1._index + 1] ||
                                (this._path._closed && t[0]))) ||
                        null
                    );
                },
                getPrevious: function () {
                    var t = this._path && this._path._curves;
                    return (
                        (t &&
                            (t[this._segment1._index - 1] ||
                                (this._path._closed && t[t.length - 1]))) ||
                        null
                    );
                },
                isSelected: function () {
                    return (
                        this.getPoint1().isSelected() &&
                        this.getHandle2().isSelected() &&
                        this.getHandle2().isSelected() &&
                        this.getPoint2().isSelected()
                    );
                },
                setSelected: function (t) {
                    this.getPoint1().setSelected(t),
                        this.getHandle1().setSelected(t),
                        this.getHandle2().setSelected(t),
                        this.getPoint2().setSelected(t);
                },
                getValues: function (t) {
                    return z.getValues(this._segment1, this._segment2, t);
                },
                getPoints: function () {
                    for (var t = this.getValues(), e = [], n = 0; 8 > n; n += 2)
                        e.push(new h(t[n], t[n + 1]));
                    return e;
                },
                getLength: function () {
                    return (
                        null == this._length &&
                            (this._length = this.isLinear()
                                ? this._segment2._point.getDistance(
                                      this._segment1._point
                                  )
                                : z.getLength(this.getValues(), 0, 1)),
                        this._length
                    );
                },
                getArea: function () {
                    return z.getArea(this.getValues());
                },
                getPart: function (t, e) {
                    return new z(z.getPart(this.getValues(), t, e));
                },
                getPartLength: function (t, e) {
                    return z.getLength(this.getValues(), t, e);
                },
                isLinear: function () {
                    return (
                        this._segment1._handleOut.isZero() &&
                        this._segment2._handleIn.isZero()
                    );
                },
                isHorizontal: function () {
                    return (
                        this.isLinear() &&
                        o.isZero(
                            this._segment1._point._y - this._segment2._point._y
                        )
                    );
                },
                getIntersections: function (t) {
                    return z.getIntersections(
                        this.getValues(),
                        t.getValues(),
                        this,
                        t,
                        []
                    );
                },
                _getParameter: function (e, n) {
                    return n
                        ? e
                        : e && e.curve === this
                        ? e.parameter
                        : e === t && n === t
                        ? 0.5
                        : this.getParameterAt(e, 0);
                },
                divide: function (t, e, n) {
                    var i = this._getParameter(t, e),
                        r = 1e-5,
                        s = null;
                    if (i > r && 1 - r > i) {
                        var a = z.subdivide(this.getValues(), i),
                            o = n ? !1 : this.isLinear(),
                            u = a[0],
                            l = a[1];
                        o ||
                            (this._segment1._handleOut.set(
                                u[2] - u[0],
                                u[3] - u[1]
                            ),
                            this._segment2._handleIn.set(
                                l[4] - l[6],
                                l[5] - l[7]
                            ));
                        var c = u[6],
                            d = u[7],
                            f = new k(
                                new h(c, d),
                                !o && new h(u[4] - c, u[5] - d),
                                !o && new h(l[2] - c, l[3] - d)
                            );
                        if (this._path)
                            this._segment1._index > 0 &&
                            0 === this._segment2._index
                                ? this._path.add(f)
                                : this._path.insert(this._segment2._index, f),
                                (s = this);
                        else {
                            var _ = this._segment2;
                            (this._segment2 = f), (s = new z(f, _));
                        }
                    }
                    return s;
                },
                split: function (t, e) {
                    return this._path
                        ? this._path.split(
                              this._segment1._index,
                              this._getParameter(t, e)
                          )
                        : null;
                },
                reverse: function () {
                    return new z(
                        this._segment2.reverse(),
                        this._segment1.reverse()
                    );
                },
                remove: function () {
                    var t = !1;
                    if (this._path) {
                        var e = this._segment2,
                            n = e._handleOut;
                        (t = e.remove()),
                            t && this._segment1._handleOut.set(n.x, n.y);
                    }
                    return t;
                },
                clone: function () {
                    return new z(this._segment1, this._segment2);
                },
                toString: function () {
                    var t = ["point1: " + this._segment1._point];
                    return (
                        this._segment1._handleOut.isZero() ||
                            t.push("handle1: " + this._segment1._handleOut),
                        this._segment2._handleIn.isZero() ||
                            t.push("handle2: " + this._segment2._handleIn),
                        t.push("point2: " + this._segment2._point),
                        "{ " + t.join(", ") + " }"
                    );
                },
                statics: {
                    getValues: function (t, e, n) {
                        var i = t._point,
                            r = t._handleOut,
                            s = e._handleIn,
                            a = e._point,
                            o = [
                                i._x,
                                i._y,
                                i._x + r._x,
                                i._y + r._y,
                                a._x + s._x,
                                a._y + s._y,
                                a._x,
                                a._y
                            ];
                        return n && n._transformCoordinates(o, o, 6), o;
                    },
                    evaluate: function (t, e, n) {
                        var i,
                            r,
                            s = t[0],
                            a = t[1],
                            o = t[2],
                            u = t[3],
                            l = t[4],
                            c = t[5],
                            d = t[6],
                            f = t[7],
                            _ = 1e-5;
                        if (0 === n && (_ > e || e > 1 - _)) {
                            var g = _ > e;
                            (i = g ? s : d), (r = g ? a : f);
                        } else {
                            var p = 3 * (o - s),
                                v = 3 * (l - o) - p,
                                m = d - s - p - v,
                                y = 3 * (u - a),
                                w = 3 * (c - u) - y,
                                x = f - a - y - w;
                            if (0 === n)
                                (i = ((m * e + v) * e + p) * e + s),
                                    (r = ((x * e + w) * e + y) * e + a);
                            else if (
                                ((_ > e && o === s && u === a) ||
                                (e > 1 - _ && l === d && c === f)
                                    ? ((i = d - s), (r = f - a))
                                    : _ > e
                                    ? ((i = p), (r = y))
                                    : e > 1 - _
                                    ? ((i = 3 * (d - l)), (r = 3 * (f - c)))
                                    : ((i = (3 * m * e + 2 * v) * e + p),
                                      (r = (3 * x * e + 2 * w) * e + y)),
                                3 === n)
                            ) {
                                var b = 6 * m * e + 2 * v,
                                    C = 6 * x * e + 2 * w;
                                return (
                                    (i * C - r * b) /
                                    Math.pow(i * i + r * r, 1.5)
                                );
                            }
                        }
                        return 2 === n ? new h(r, -i) : new h(i, r);
                    },
                    subdivide: function (e, n) {
                        var i = e[0],
                            r = e[1],
                            s = e[2],
                            a = e[3],
                            o = e[4],
                            h = e[5],
                            u = e[6],
                            l = e[7];
                        n === t && (n = 0.5);
                        var c = 1 - n,
                            d = c * i + n * s,
                            f = c * r + n * a,
                            _ = c * s + n * o,
                            g = c * a + n * h,
                            p = c * o + n * u,
                            v = c * h + n * l,
                            m = c * d + n * _,
                            y = c * f + n * g,
                            w = c * _ + n * p,
                            x = c * g + n * v,
                            b = c * m + n * w,
                            C = c * y + n * x;
                        return [
                            [i, r, d, f, m, y, b, C],
                            [b, C, w, x, p, v, u, l]
                        ];
                    },
                    solveCubic: function (t, e, n, i, r, s) {
                        var a = t[e],
                            h = t[e + 2],
                            u = t[e + 4],
                            l = t[e + 6],
                            c = 3 * (h - a),
                            d = 3 * (u - h) - c,
                            f = l - a - c - d;
                        return o.solveCubic(f, d, c, a - n, i, r, s);
                    },
                    getParameterOf: function (t, e, n) {
                        var i = 1e-5;
                        if (Math.abs(t[0] - e) < i && Math.abs(t[1] - n) < i)
                            return 0;
                        if (Math.abs(t[6] - e) < i && Math.abs(t[7] - n) < i)
                            return 1;
                        for (
                            var r,
                                s,
                                a = [],
                                o = [],
                                h = z.solveCubic(t, 0, e, a),
                                u = z.solveCubic(t, 1, n, o),
                                l = 0;
                            -1 == h || h > l;

                        )
                            if (-1 == h || ((r = a[l++]) >= 0 && 1 >= r)) {
                                for (var c = 0; -1 == u || u > c; )
                                    if (
                                        (-1 == u ||
                                            ((s = o[c++]) >= 0 && 1 >= s)) &&
                                        (-1 == h ? (r = s) : -1 == u && (s = r),
                                        Math.abs(r - s) < i)
                                    )
                                        return 0.5 * (r + s);
                                if (-1 == h) break;
                            }
                        return null;
                    },
                    getPart: function (t, e, n) {
                        return (
                            e > 0 && (t = z.subdivide(t, e)[1]),
                            1 > n && (t = z.subdivide(t, (n - e) / (1 - e))[0]),
                            t
                        );
                    },
                    isLinear: function (t) {
                        var e = o.isZero;
                        return (
                            e(t[0] - t[2]) &&
                            e(t[1] - t[3]) &&
                            e(t[4] - t[6]) &&
                            e(t[5] - t[7])
                        );
                    },
                    isFlatEnough: function (t, e) {
                        var n = t[0],
                            i = t[1],
                            r = t[2],
                            s = t[3],
                            a = t[4],
                            o = t[5],
                            h = t[6],
                            u = t[7],
                            l = 3 * r - 2 * n - h,
                            c = 3 * s - 2 * i - u,
                            d = 3 * a - 2 * h - n,
                            f = 3 * o - 2 * u - i;
                        return (
                            Math.max(l * l, d * d) + Math.max(c * c, f * f) <
                            10 * e * e
                        );
                    },
                    getArea: function (t) {
                        var e = t[0],
                            n = t[1],
                            i = t[2],
                            r = t[3],
                            s = t[4],
                            a = t[5],
                            o = t[6],
                            h = t[7];
                        return (
                            (3 * r * e -
                                1.5 * r * s -
                                1.5 * r * o -
                                3 * n * i -
                                1.5 * n * s -
                                0.5 * n * o +
                                1.5 * a * e +
                                1.5 * a * i -
                                3 * a * o +
                                0.5 * h * e +
                                1.5 * h * i +
                                3 * h * s) /
                            10
                        );
                    },
                    getBounds: function (t) {
                        for (
                            var e = t.slice(0, 2),
                                n = e.slice(),
                                i = [0, 0],
                                r = 0;
                            2 > r;
                            r++
                        )
                            z._addBounds(
                                t[r],
                                t[r + 2],
                                t[r + 4],
                                t[r + 6],
                                r,
                                0,
                                e,
                                n,
                                i
                            );
                        return new f(e[0], e[1], n[0] - e[0], n[1] - e[1]);
                    },
                    _addBounds: function (t, e, n, i, r, s, a, h, u) {
                        function l(t, e) {
                            var n = t - e,
                                i = t + e;
                            n < a[r] && (a[r] = n), i > h[r] && (h[r] = i);
                        }
                        var c = 3 * (e - n) - t + i,
                            d = 2 * (t + n) - 4 * e,
                            f = e - t,
                            _ = o.solveQuadratic(c, d, f, u),
                            g = 1e-5,
                            p = 1 - g;
                        l(i, 0);
                        for (var v = 0; _ > v; v++) {
                            var m = u[v],
                                y = 1 - m;
                            m > g &&
                                p > m &&
                                l(
                                    y * y * y * t +
                                        3 * y * y * m * e +
                                        3 * y * m * m * n +
                                        m * m * m * i,
                                    s
                                );
                        }
                    }
                }
            },
            e.each(
                [
                    "getBounds",
                    "getStrokeBounds",
                    "getHandleBounds",
                    "getRoughBounds"
                ],
                function (t) {
                    this[t] = function () {
                        this._bounds || (this._bounds = {});
                        var e = this._bounds[t];
                        return (
                            e ||
                                (e = this._bounds[t] =
                                    O[t](
                                        [this._segment1, this._segment2],
                                        !1,
                                        this._path.getStyle()
                                    )),
                            e.clone()
                        );
                    };
                },
                {}
            ),
            e.each(
                ["getPoint", "getTangent", "getNormal", "getCurvature"],
                function (t, e) {
                    (this[t + "At"] = function (t, n) {
                        var i = this.getValues();
                        return z.evaluate(
                            i,
                            n ? t : z.getParameterAt(i, t, 0),
                            e
                        );
                    }),
                        (this[t] = function (t) {
                            return z.evaluate(this.getValues(), t, e);
                        });
                },
                {
                    beans: !1,
                    getParameterAt: function (e, n) {
                        return z.getParameterAt(
                            this.getValues(),
                            e,
                            n !== t ? n : 0 > e ? 1 : 0
                        );
                    },
                    getParameterOf: function () {
                        var t = h.read(arguments);
                        return z.getParameterOf(this.getValues(), t.x, t.y);
                    },
                    getLocationAt: function (t, e) {
                        return (
                            e || (t = this.getParameterAt(t)), new A(this, t)
                        );
                    },
                    getLocationOf: function () {
                        var t = h.read(arguments),
                            e = this.getParameterOf(t);
                        return null != e ? new A(this, e) : null;
                    },
                    getOffsetOf: function () {
                        var t = this.getLocationOf.apply(this, arguments);
                        return t ? t.getOffset() : null;
                    },
                    getNearestLocation: function () {
                        function t(t) {
                            if (t >= 0 && 1 >= t) {
                                var i = e.getDistance(z.evaluate(n, t, 0), !0);
                                if (r > i) return (r = i), (s = t), !0;
                            }
                        }
                        for (
                            var e = h.read(arguments),
                                n = this.getValues(),
                                i = 100,
                                r = 1 / 0,
                                s = 0,
                                a = 0;
                            i >= a;
                            a++
                        )
                            t(a / i);
                        for (var o = 1 / (2 * i); o > 1e-5; )
                            t(s - o) || t(s + o) || (o /= 2);
                        var u = z.evaluate(n, s, 0);
                        return new A(
                            this,
                            s,
                            u,
                            null,
                            null,
                            null,
                            e.getDistance(u)
                        );
                    },
                    getNearestPoint: function () {
                        return this.getNearestLocation
                            .apply(this, arguments)
                            .getPoint();
                    }
                }
            ),
            new (function () {
                function e(t) {
                    var e = t[0],
                        n = t[1],
                        i = t[2],
                        r = t[3],
                        s = t[4],
                        a = t[5],
                        o = t[6],
                        h = t[7],
                        u = 9 * (i - s) + 3 * (o - e),
                        l = 6 * (e + s) - 12 * i,
                        c = 3 * (i - e),
                        d = 9 * (r - a) + 3 * (h - n),
                        f = 6 * (n + a) - 12 * r,
                        _ = 3 * (r - n);
                    return function (t) {
                        var e = (u * t + l) * t + c,
                            n = (d * t + f) * t + _;
                        return Math.sqrt(e * e + n * n);
                    };
                }
                function n(t, e) {
                    return Math.max(
                        2,
                        Math.min(16, Math.ceil(32 * Math.abs(e - t)))
                    );
                }
                return {
                    statics: !0,
                    getLength: function (i, r, s) {
                        r === t && (r = 0), s === t && (s = 1);
                        var a = o.isZero;
                        if (
                            0 === r &&
                            1 === s &&
                            a(i[0] - i[2]) &&
                            a(i[1] - i[3]) &&
                            a(i[6] - i[4]) &&
                            a(i[7] - i[5])
                        ) {
                            var h = i[6] - i[0],
                                u = i[7] - i[1];
                            return Math.sqrt(h * h + u * u);
                        }
                        var l = e(i);
                        return o.integrate(l, r, s, n(r, s));
                    },
                    getParameterAt: function (t, i, r) {
                        function s(t) {
                            var e = n(r, t);
                            return (
                                (f +=
                                    t > r
                                        ? o.integrate(l, r, t, e)
                                        : -o.integrate(l, t, r, e)),
                                (r = t),
                                f - i
                            );
                        }
                        if (0 === i) return r;
                        var a = i > 0,
                            h = a ? r : 0,
                            u = a ? 1 : r,
                            i = Math.abs(i),
                            l = e(t),
                            c = o.integrate(l, h, u, n(h, u));
                        if (i >= c) return a ? u : h;
                        var d = i / c,
                            f = 0;
                        return o.findRoot(
                            s,
                            l,
                            a ? h + d : u - d,
                            h,
                            u,
                            16,
                            1e-5
                        );
                    }
                };
            })(),
            new (function () {
                function t(t, e, n, i, r, s, a, o) {
                    var h = new A(n, i, r, s, a, o);
                    (!e || e(h)) && t.push(h);
                }
                function e(r, s, a, o, h, u, l, c, d, f, _, g, v) {
                    if (!(v > 20)) {
                        var m,
                            y,
                            w,
                            x = s[0],
                            b = s[1],
                            C = s[6],
                            S = s[7],
                            P = 1e-5,
                            k = 1e-9,
                            M = p.getSignedDistance,
                            A = M(x, b, C, S, s[2], s[3]) || 0,
                            I = M(x, b, C, S, s[4], s[5]) || 0,
                            O = A * I > 0 ? 0.75 : 4 / 9,
                            T = O * Math.min(0, A, I),
                            L = O * Math.max(0, A, I),
                            E = M(x, b, C, S, r[0], r[1]),
                            N = M(x, b, C, S, r[2], r[3]),
                            j = M(x, b, C, S, r[4], r[5]),
                            D = M(x, b, C, S, r[6], r[7]);
                        if (x === C && k >= f - d && v > 3)
                            (m = (c + l) / 2), (y = m), (w = 0);
                        else {
                            var B,
                                R,
                                F = n(E, N, j, D),
                                V = F[0],
                                q = F[1];
                            if (
                                ((B = i(V, q, T, L)),
                                V.reverse(),
                                q.reverse(),
                                (R = i(V, q, T, L)),
                                null == B || null == R)
                            )
                                return !1;
                            (r = z.getPart(r, B, R)),
                                (w = R - B),
                                (m = c * B + l * (1 - B)),
                                (y = c * R + l * (1 - R));
                        }
                        if (_ > 0.8 && w > 0.8)
                            if (y - m > f - d) {
                                var Z = z.subdivide(r, 0.5),
                                    U = m + (y - m) / 2;
                                e(s, Z[0], o, a, h, u, d, f, m, U, w, !g, ++v),
                                    e(
                                        s,
                                        Z[1],
                                        o,
                                        a,
                                        h,
                                        u,
                                        d,
                                        f,
                                        U,
                                        y,
                                        w,
                                        !g,
                                        v
                                    );
                            } else {
                                var Z = z.subdivide(s, 0.5),
                                    U = d + (f - d) / 2;
                                e(Z[0], r, o, a, h, u, d, U, m, y, w, !g, ++v),
                                    e(
                                        Z[1],
                                        r,
                                        o,
                                        a,
                                        h,
                                        u,
                                        U,
                                        f,
                                        m,
                                        y,
                                        w,
                                        !g,
                                        v
                                    );
                            }
                        else if (Math.max(f - d, y - m) < P) {
                            var H = m + (y - m) / 2,
                                W = d + (f - d) / 2;
                            g
                                ? t(
                                      h,
                                      u,
                                      o,
                                      W,
                                      z.evaluate(s, W, 0),
                                      a,
                                      H,
                                      z.evaluate(r, H, 0)
                                  )
                                : t(
                                      h,
                                      u,
                                      a,
                                      H,
                                      z.evaluate(r, H, 0),
                                      o,
                                      W,
                                      z.evaluate(s, W, 0)
                                  );
                        } else e(s, r, o, a, h, u, d, f, m, y, w, !g, ++v);
                    }
                }
                function n(t, e, n, i) {
                    var r,
                        s = [0, t],
                        a = [1 / 3, e],
                        o = [2 / 3, n],
                        h = [1, i],
                        u = p.getSignedDistance,
                        l = u(0, t, 1, i, 1 / 3, e),
                        c = u(0, t, 1, i, 2 / 3, n),
                        d = !1;
                    if (0 > l * c)
                        (r = [
                            [s, a, h],
                            [s, o, h]
                        ]),
                            (d = 0 > l);
                    else {
                        var f,
                            _ = 0,
                            g = 0 === l || 0 === c;
                        Math.abs(l) > Math.abs(c)
                            ? ((f = a),
                              (_ =
                                  ((i - n - (i - t) / 3) *
                                      (2 * (i - n) - i + e)) /
                                  3))
                            : ((f = o),
                              (_ =
                                  ((e - t + (t - i) / 3) *
                                      (-2 * (t - e) + t - n)) /
                                  3)),
                            (r =
                                0 > _ || g
                                    ? [
                                          [s, f, h],
                                          [s, h]
                                      ]
                                    : [
                                          [s, a, o, h],
                                          [s, h]
                                      ]),
                            (d = l ? 0 > l : 0 > c);
                    }
                    return d ? r.reverse() : r;
                }
                function i(t, e, n, i) {
                    for (
                        var r, s, a, o, h, u = null, l = 0, c = e.length - 1;
                        c > l;
                        l++
                    ) {
                        if (((a = e[l][1]), (h = e[l + 1][1]), h > a)) r = null;
                        else {
                            if (!(i >= h)) continue;
                            (s = e[l][0]),
                                (o = e[l + 1][0]),
                                (r = s + ((i - a) * (o - s)) / (h - a));
                        }
                        break;
                    }
                    t[0][1] <= i && (r = t[0][0]);
                    for (var l = 0, c = t.length - 1; c > l; l++) {
                        if (((a = t[l][1]), (h = t[l + 1][1]), a >= n)) u = r;
                        else if (a > h) u = null;
                        else {
                            if (!(h >= n)) continue;
                            (s = t[l][0]),
                                (o = t[l + 1][0]),
                                (u = s + ((n - a) * (o - s)) / (h - a));
                        }
                        break;
                    }
                    return u;
                }
                function r(e, n, i, r, s, a) {
                    for (
                        var o = z.isLinear(e),
                            h = o ? n : e,
                            u = o ? e : n,
                            l = u[0],
                            c = u[1],
                            d = u[6],
                            f = u[7],
                            _ = d - l,
                            g = f - c,
                            p = Math.atan2(-g, _),
                            v = Math.sin(p),
                            m = Math.cos(p),
                            y = _ * m - g * v,
                            w = [0, 0, 0, 0, y, 0, y, 0],
                            x = [],
                            b = 0;
                        8 > b;
                        b += 2
                    ) {
                        var C = h[b] - l,
                            S = h[b + 1] - c;
                        x.push(C * m - S * v, S * m + C * v);
                    }
                    for (
                        var P = [], k = z.solveCubic(x, 1, 0, P, 0, 1), b = 0;
                        k > b;
                        b++
                    ) {
                        var M = P[b],
                            C = z.evaluate(x, M, 0).x;
                        if (C >= 0 && y >= C) {
                            var A = z.getParameterOf(w, C, 0),
                                I = o ? A : M,
                                O = o ? M : A;
                            t(
                                s,
                                a,
                                i,
                                I,
                                z.evaluate(e, I, 0),
                                r,
                                O,
                                z.evaluate(n, O, 0)
                            );
                        }
                    }
                }
                function s(e, n, i, r, s, a) {
                    var o = p.intersect(
                        e[0],
                        e[1],
                        e[6],
                        e[7],
                        n[0],
                        n[1],
                        n[6],
                        n[7]
                    );
                    if (o) {
                        var h = o.x,
                            u = o.y;
                        t(
                            s,
                            a,
                            i,
                            z.getParameterOf(e, h, u),
                            o,
                            r,
                            z.getParameterOf(n, h, u),
                            o
                        );
                    }
                }
                return {
                    statics: {
                        getIntersections: function (t, n, i, a, o, h) {
                            var u = z.isLinear(t),
                                l = z.isLinear(n);
                            return (
                                (u && l ? s : u || l ? r : e)(
                                    t,
                                    n,
                                    i,
                                    a,
                                    o,
                                    h,
                                    0,
                                    1,
                                    0,
                                    1,
                                    0,
                                    !1,
                                    0
                                ),
                                o
                            );
                        }
                    }
                };
            })()
        ),
        A = e.extend(
            {
                _class: "CurveLocation",
                beans: !0,
                initialize: function me(t, e, n, i, r, s, a) {
                    (this._id = me._id = (me._id || 0) + 1),
                        (this._curve = t),
                        (this._segment1 = t._segment1),
                        (this._segment2 = t._segment2),
                        (this._parameter = e),
                        (this._point = n),
                        (this._curve2 = i),
                        (this._parameter2 = r),
                        (this._point2 = s),
                        (this._distance = a);
                },
                getSegment: function (t) {
                    if (!this._segment) {
                        var e = this.getCurve(),
                            n = this.getParameter();
                        if (1 === n) this._segment = e._segment2;
                        else if (0 === n || t) this._segment = e._segment1;
                        else {
                            if (null == n) return null;
                            this._segment =
                                e.getPartLength(0, n) < e.getPartLength(n, 1)
                                    ? e._segment1
                                    : e._segment2;
                        }
                    }
                    return this._segment;
                },
                getCurve: function (t) {
                    return (
                        (!this._curve || t) &&
                            ((this._curve = this._segment1.getCurve()),
                            null == this._curve.getParameterOf(this._point) &&
                                (this._curve = this._segment2
                                    .getPrevious()
                                    .getCurve())),
                        this._curve
                    );
                },
                getIntersection: function () {
                    var t = this._intersection;
                    if (!t && this._curve2) {
                        var e = this._parameter2;
                        (this._intersection = t =
                            new A(
                                this._curve2,
                                e,
                                this._point2 || this._point,
                                this
                            )),
                            (t._intersection = this);
                    }
                    return t;
                },
                getPath: function () {
                    var t = this.getCurve();
                    return t && t._path;
                },
                getIndex: function () {
                    var t = this.getCurve();
                    return t && t.getIndex();
                },
                getOffset: function () {
                    var t = this.getPath();
                    return t ? t._getOffset(this) : this.getCurveOffset();
                },
                getCurveOffset: function () {
                    var t = this.getCurve(),
                        e = this.getParameter();
                    return null != e && t && t.getPartLength(0, e);
                },
                getParameter: function (t) {
                    if ((null == this._parameter || t) && this._point) {
                        var e = this.getCurve(t && this._point);
                        this._parameter = e && e.getParameterOf(this._point);
                    }
                    return this._parameter;
                },
                getPoint: function (t) {
                    if ((!this._point || t) && null != this._parameter) {
                        var e = this.getCurve();
                        this._point = e && e.getPointAt(this._parameter, !0);
                    }
                    return this._point;
                },
                getDistance: function () {
                    return this._distance;
                },
                divide: function () {
                    var t = this.getCurve(!0);
                    return t && t.divide(this.getParameter(!0), !0);
                },
                split: function () {
                    var t = this.getCurve(!0);
                    return t && t.split(this.getParameter(!0), !0);
                },
                equals: function (t) {
                    var e = o.isZero;
                    return (
                        this === t ||
                        (t &&
                            this._curve === t._curve &&
                            this._curve2 === t._curve2 &&
                            e(this._parameter - t._parameter) &&
                            e(this._parameter2 - t._parameter2)) ||
                        !1
                    );
                },
                toString: function () {
                    var t = [],
                        e = this.getPoint(),
                        n = a.instance;
                    e && t.push("point: " + e);
                    var i = this.getIndex();
                    null != i && t.push("index: " + i);
                    var r = this.getParameter();
                    return (
                        null != r && t.push("parameter: " + n.number(r)),
                        null != this._distance &&
                            t.push("distance: " + n.number(this._distance)),
                        "{ " + t.join(", ") + " }"
                    );
                }
            },
            e.each(
                ["Tangent", "Normal", "Curvature"],
                function (t) {
                    var e = "get" + t + "At";
                    this["get" + t] = function () {
                        var t = this.getParameter(),
                            n = this.getCurve();
                        return null != t && n && n[e](t, !0);
                    };
                },
                {}
            )
        ),
        I = y.extend({
            _class: "PathItem",
            initialize: function () {},
            getIntersections: function (e, n) {
                function i(t, e) {
                    var n = t.getPath(),
                        i = e.getPath();
                    return n === i
                        ? t.getIndex() +
                              t.getParameter() -
                              (e.getIndex() + e.getParameter())
                        : n._id - i._id;
                }
                if (
                    (this === e && (e = null),
                    e && !this.getBounds().touches(e.getBounds()))
                )
                    return [];
                for (
                    var r = [],
                        s = this.getCurves(),
                        a = e ? e.getCurves() : s,
                        o = this._matrix.orNullIfIdentity(),
                        h = e ? e._matrix.orNullIfIdentity() : o,
                        u = s.length,
                        l = e ? a.length : u,
                        c = [],
                        d = 1e-11,
                        f = 1 - 1e-11,
                        _ = 0;
                    l > _;
                    _++
                )
                    c[_] = a[_].getValues(h);
                for (var _ = 0; u > _; _++) {
                    var g = s[_],
                        v = e ? g.getValues(o) : c[_];
                    if (!e) {
                        var m = g.getSegment1(),
                            y = g.getSegment2(),
                            w = m._handleOut,
                            x = y._handleIn;
                        if (
                            new p(
                                m._point.subtract(w),
                                w.multiply(2),
                                !0
                            ).intersect(
                                new p(y._point.subtract(x), x.multiply(2), !0),
                                !1
                            )
                        ) {
                            var b = z.subdivide(v);
                            z.getIntersections(
                                b[0],
                                b[1],
                                g,
                                g,
                                r,
                                function (e) {
                                    return e._parameter <= f
                                        ? ((e._parameter /= 2),
                                          (e._parameter2 =
                                              0.5 + e._parameter2 / 2),
                                          !0)
                                        : t;
                                }
                            );
                        }
                    }
                    for (var C = e ? 0 : _ + 1; l > C; C++)
                        z.getIntersections(
                            v,
                            c[C],
                            g,
                            a[C],
                            r,
                            !e &&
                                (C === _ + 1 || (C === l - 1 && 0 === _)) &&
                                function (t) {
                                    var e = t._parameter;
                                    return e >= d && f >= e;
                                }
                        );
                }
                for (var S = r.length - 1, _ = S; _ >= 0; _--) {
                    var P = r[_],
                        k = P._curve.getNext(),
                        M = P._curve2.getNext();
                    k &&
                        P._parameter >= f &&
                        ((P._parameter = 0), (P._curve = k)),
                        M &&
                            P._parameter2 >= f &&
                            ((P._parameter2 = 0), (P._curve2 = M));
                }
                if (S > 0) {
                    r.sort(i);
                    for (var _ = S; _ >= 1; _--)
                        r[_].equals(r[0 === _ ? S : _ - 1]) &&
                            (r.splice(_, 1), S--);
                }
                if (n) {
                    for (var _ = S; _ >= 0; _--) r.push(r[_].getIntersection());
                    r.sort(i);
                }
                return r;
            },
            setPathData: function (t) {
                function e(t, e) {
                    var n = +i[t];
                    return o && (n += u[e]), n;
                }
                function n(t) {
                    return new h(e(t, "x"), e(t + 1, "y"));
                }
                var i,
                    r,
                    s,
                    a = t.match(/[mlhvcsqtaz][^mlhvcsqtaz]*/gi),
                    o = !1,
                    u = new h(),
                    l = new h();
                this.clear();
                for (var d = 0, f = a.length; f > d; d++) {
                    var _ = a[d],
                        g = _[0],
                        p = g.toLowerCase();
                    i = _.match(/[+-]?(?:\d*\.\d+|\d+\.?)(?:[eE][+-]?\d+)?/g);
                    var v = i && i.length;
                    switch (
                        ((o = g === p),
                        "z" !== r || /[mz]/.test(p) || this.moveTo((u = l)),
                        p)
                    ) {
                        case "m":
                        case "l":
                            var m = "m" === p;
                            m && r && "z" !== r && this.closePath(!0);
                            for (var y = 0; v > y; y += 2)
                                this[0 === y && m ? "moveTo" : "lineTo"](
                                    (u = n(y))
                                );
                            (s = u), m && (l = u);
                            break;
                        case "h":
                        case "v":
                            for (
                                var w = "h" === p ? "x" : "y", y = 0;
                                v > y;
                                y++
                            )
                                (u[w] = e(y, w)), this.lineTo(u);
                            s = u;
                            break;
                        case "c":
                            for (var y = 0; v > y; y += 6)
                                this.cubicCurveTo(
                                    n(y),
                                    (s = n(y + 2)),
                                    (u = n(y + 4))
                                );
                            break;
                        case "s":
                            for (var y = 0; v > y; y += 4)
                                this.cubicCurveTo(
                                    /[cs]/.test(r)
                                        ? u.multiply(2).subtract(s)
                                        : u,
                                    (s = n(y)),
                                    (u = n(y + 2))
                                ),
                                    (r = p);
                            break;
                        case "q":
                            for (var y = 0; v > y; y += 4)
                                this.quadraticCurveTo(
                                    (s = n(y)),
                                    (u = n(y + 2))
                                );
                            break;
                        case "t":
                            for (var y = 0; v > y; y += 2)
                                this.quadraticCurveTo(
                                    (s = /[qt]/.test(r)
                                        ? u.multiply(2).subtract(s)
                                        : u),
                                    (u = n(y))
                                ),
                                    (r = p);
                            break;
                        case "a":
                            for (var y = 0; v > y; y += 7)
                                this.arcTo(
                                    (u = n(y + 5)),
                                    new c(+i[0], +i[1]),
                                    +i[2],
                                    +i[4],
                                    +i[3]
                                );
                            break;
                        case "z":
                            this.closePath(!0);
                    }
                    r = p;
                }
            },
            _canComposite: function () {
                return !(this.hasFill() && this.hasStroke());
            },
            _contains: function (t) {
                var e = this._getWinding(t, !1, !0);
                return !!("evenodd" === this.getWindingRule() ? 1 & e : e);
            }
        }),
        O = I.extend(
            {
                _class: "Path",
                _serializeFields: {
                    segments: [],
                    closed: !1
                },
                initialize: function (e) {
                    (this._closed = !1), (this._segments = []);
                    var n = Array.isArray(e)
                        ? "object" == typeof e[0]
                            ? e
                            : arguments
                        : !e || e.size !== t || (e.x === t && e.point === t)
                        ? null
                        : arguments;
                    n && n.length > 0
                        ? this.setSegments(n)
                        : ((this._curves = t),
                          (this._selectedSegmentState = 0),
                          n ||
                              "string" != typeof e ||
                              (this.setPathData(e), (e = null))),
                        this._initialize(!n && e);
                },
                _equals: function (t) {
                    return e.equals(this._segments, t._segments);
                },
                clone: function (e) {
                    var n = new O(y.NO_INSERT);
                    return (
                        n.setSegments(this._segments),
                        (n._closed = this._closed),
                        this._clockwise !== t &&
                            (n._clockwise = this._clockwise),
                        this._clone(n, e)
                    );
                },
                _changed: function ye(e) {
                    if ((ye.base.call(this, e), 8 & e)) {
                        var n = this._parent;
                        if (
                            (n && (n._currentPath = t),
                            (this._length = this._clockwise = t),
                            this._curves && !(16 & e))
                        )
                            for (var i = 0, r = this._curves.length; r > i; i++)
                                this._curves[i]._changed();
                        this._monoCurves = t;
                    } else 32 & e && (this._bounds = t);
                },
                getStyle: function () {
                    var t = this._parent;
                    return (t instanceof T ? t : this)._style;
                },
                getSegments: function () {
                    return this._segments;
                },
                setSegments: function (e) {
                    var n = this.isFullySelected();
                    (this._segments.length = 0),
                        (this._selectedSegmentState = 0),
                        (this._curves = t),
                        e && e.length > 0 && this._add(k.readAll(e)),
                        n && this.setFullySelected(!0);
                },
                getFirstSegment: function () {
                    return this._segments[0];
                },
                getLastSegment: function () {
                    return this._segments[this._segments.length - 1];
                },
                getCurves: function () {
                    var t = this._curves,
                        e = this._segments;
                    if (!t) {
                        var n = this._countCurves();
                        t = this._curves = Array(n);
                        for (var i = 0; n > i; i++)
                            t[i] = new z(this, e[i], e[i + 1] || e[0]);
                    }
                    return t;
                },
                getFirstCurve: function () {
                    return this.getCurves()[0];
                },
                getLastCurve: function () {
                    var t = this.getCurves();
                    return t[t.length - 1];
                },
                isClosed: function () {
                    return this._closed;
                },
                setClosed: function (t) {
                    if (this._closed != (t = !!t)) {
                        if (((this._closed = t), this._curves)) {
                            var e = (this._curves.length = this._countCurves());
                            t &&
                                (this._curves[e - 1] = new z(
                                    this,
                                    this._segments[e - 1],
                                    this._segments[0]
                                ));
                        }
                        this._changed(25);
                    }
                }
            },
            {
                beans: !0,
                getPathData: function (t, e) {
                    function n(e, n) {
                        e._transformCoordinates(t, g, !1),
                            (i = g[0]),
                            (r = g[1]),
                            p
                                ? (v.push("M" + _.pair(i, r)), (p = !1))
                                : ((h = g[2]),
                                  (u = g[3]),
                                  h === i && u === r && l === s && c === o
                                      ? n || v.push("l" + _.pair(i - s, r - o))
                                      : v.push(
                                            "c" +
                                                _.pair(l - s, c - o) +
                                                " " +
                                                _.pair(h - s, u - o) +
                                                " " +
                                                _.pair(i - s, r - o)
                                        )),
                            (s = i),
                            (o = r),
                            (l = g[4]),
                            (c = g[5]);
                    }
                    var i,
                        r,
                        s,
                        o,
                        h,
                        u,
                        l,
                        c,
                        d = this._segments,
                        f = d.length,
                        _ = new a(e),
                        g = Array(6),
                        p = !0,
                        v = [];
                    if (0 === f) return "";
                    for (var m = 0; f > m; m++) n(d[m]);
                    return (
                        this._closed && f > 0 && (n(d[0], !0), v.push("z")),
                        v.join("")
                    );
                }
            },
            {
                isEmpty: function () {
                    return 0 === this._segments.length;
                },
                isPolygon: function () {
                    for (var t = 0, e = this._segments.length; e > t; t++)
                        if (!this._segments[t].isLinear()) return !1;
                    return !0;
                },
                _transformContent: function (t) {
                    for (
                        var e = Array(6), n = 0, i = this._segments.length;
                        i > n;
                        n++
                    )
                        this._segments[n]._transformCoordinates(t, e, !0);
                    return !0;
                },
                _add: function (t, e) {
                    for (
                        var n = this._segments,
                            i = this._curves,
                            r = t.length,
                            s = null == e,
                            e = s ? n.length : e,
                            a = 0;
                        r > a;
                        a++
                    ) {
                        var o = t[a];
                        o._path && (o = t[a] = o.clone()),
                            (o._path = this),
                            (o._index = e + a),
                            o._selectionState &&
                                this._updateSelection(o, 0, o._selectionState);
                    }
                    if (s) n.push.apply(n, t);
                    else {
                        n.splice.apply(n, [e, 0].concat(t));
                        for (var a = e + r, h = n.length; h > a; a++)
                            n[a]._index = a;
                    }
                    if (i || t._curves) {
                        i || (i = this._curves = []);
                        var u = e > 0 ? e - 1 : e,
                            l = u,
                            c = Math.min(u + r, this._countCurves());
                        t._curves &&
                            (i.splice.apply(i, [u, 0].concat(t._curves)),
                            (l += t._curves.length));
                        for (var a = l; c > a; a++)
                            i.splice(a, 0, new z(this, null, null));
                        this._adjustCurves(u, c);
                    }
                    return this._changed(25), t;
                },
                _adjustCurves: function (t, e) {
                    for (
                        var n, i = this._segments, r = this._curves, s = t;
                        e > s;
                        s++
                    )
                        (n = r[s]),
                            (n._path = this),
                            (n._segment1 = i[s]),
                            (n._segment2 = i[s + 1] || i[0]),
                            n._changed();
                    (n = r[this._closed && 0 === t ? i.length - 1 : t - 1]) &&
                        ((n._segment2 = i[t] || i[0]), n._changed()),
                        (n = r[e]) && ((n._segment1 = i[e]), n._changed());
                },
                _countCurves: function () {
                    var t = this._segments.length;
                    return !this._closed && t > 0 ? t - 1 : t;
                },
                add: function (t) {
                    return arguments.length > 1 && "number" != typeof t
                        ? this._add(k.readAll(arguments))
                        : this._add([k.read(arguments)])[0];
                },
                insert: function (t, e) {
                    return arguments.length > 2 && "number" != typeof e
                        ? this._add(k.readAll(arguments, 1), t)
                        : this._add([k.read(arguments, 1)], t)[0];
                },
                addSegment: function () {
                    return this._add([k.read(arguments)])[0];
                },
                insertSegment: function (t) {
                    return this._add([k.read(arguments, 1)], t)[0];
                },
                addSegments: function (t) {
                    return this._add(k.readAll(t));
                },
                insertSegments: function (t, e) {
                    return this._add(k.readAll(e), t);
                },
                removeSegment: function (t) {
                    return this.removeSegments(t, t + 1)[0] || null;
                },
                removeSegments: function (t, n, i) {
                    (t = t || 0), (n = e.pick(n, this._segments.length));
                    var r = this._segments,
                        s = this._curves,
                        a = r.length,
                        o = r.splice(t, n - t),
                        h = o.length;
                    if (!h) return o;
                    for (var u = 0; h > u; u++) {
                        var l = o[u];
                        l._selectionState &&
                            this._updateSelection(l, l._selectionState, 0),
                            (l._index = l._path = null);
                    }
                    for (var u = t, c = r.length; c > u; u++) r[u]._index = u;
                    if (s) {
                        var d =
                                t > 0 && n === a + (this._closed ? 1 : 0)
                                    ? t - 1
                                    : t,
                            s = s.splice(d, h);
                        i && (o._curves = s.slice(1)), this._adjustCurves(d, d);
                    }
                    return this._changed(25), o;
                },
                clear: "#removeSegments",
                getLength: function () {
                    if (null == this._length) {
                        var t = this.getCurves();
                        this._length = 0;
                        for (var e = 0, n = t.length; n > e; e++)
                            this._length += t[e].getLength();
                    }
                    return this._length;
                },
                getArea: function () {
                    for (
                        var t = this.getCurves(), e = 0, n = 0, i = t.length;
                        i > n;
                        n++
                    )
                        e += t[n].getArea();
                    return e;
                },
                isFullySelected: function () {
                    var t = this._segments.length;
                    return (
                        this._selected &&
                        t > 0 &&
                        this._selectedSegmentState === 7 * t
                    );
                },
                setFullySelected: function (t) {
                    t && this._selectSegments(!0), this.setSelected(t);
                },
                setSelected: function we(t) {
                    t || this._selectSegments(!1), we.base.call(this, t);
                },
                _selectSegments: function (t) {
                    var e = this._segments.length;
                    this._selectedSegmentState = t ? 7 * e : 0;
                    for (var n = 0; e > n; n++)
                        this._segments[n]._selectionState = t ? 7 : 0;
                },
                _updateSelection: function (t, e, n) {
                    t._selectionState = n;
                    var i = (this._selectedSegmentState += n - e);
                    i > 0 && this.setSelected(!0);
                },
                flatten: function (t) {
                    for (
                        var e = new L(this),
                            n = 0,
                            i = e.length / Math.ceil(e.length / t),
                            r = e.length + (this._closed ? -i : i) / 2,
                            s = [];
                        r >= n;

                    )
                        s.push(new k(e.evaluate(n, 0))), (n += i);
                    this.setSegments(s);
                },
                reduce: function () {
                    for (
                        var t = this.getCurves(), e = t.length - 1;
                        e >= 0;
                        e--
                    ) {
                        var n = t[e];
                        n.isLinear() && 0 === n.getLength() && n.remove();
                    }
                    return this;
                },
                simplify: function (t) {
                    if (this._segments.length > 2) {
                        var e = new E(this, t || 2.5);
                        this.setSegments(e.fit());
                    }
                },
                split: function (t, e) {
                    if (null !== e) {
                        if (1 === arguments.length) {
                            var n = t;
                            "number" == typeof n && (n = this.getLocationAt(n)),
                                (t = n.index),
                                (e = n.parameter);
                        }
                        var i = 1e-5;
                        e >= 1 - i && (t++, e--);
                        var r = this.getCurves();
                        if (t >= 0 && t < r.length) {
                            e > i && r[t++].divide(e, !0);
                            var s,
                                a = this.removeSegments(
                                    t,
                                    this._segments.length,
                                    !0
                                );
                            return (
                                this._closed
                                    ? (this.setClosed(!1), (s = this))
                                    : t > 0 &&
                                      (s = this._clone(
                                          new O().insertAbove(this, !0)
                                      )),
                                s._add(a, 0),
                                this.addSegment(a[0]),
                                s
                            );
                        }
                        return null;
                    }
                },
                isClockwise: function () {
                    return this._clockwise !== t
                        ? this._clockwise
                        : O.isClockwise(this._segments);
                },
                setClockwise: function (t) {
                    this.isClockwise() != (t = !!t) && this.reverse(),
                        (this._clockwise = t);
                },
                reverse: function () {
                    this._segments.reverse();
                    for (var e = 0, n = this._segments.length; n > e; e++) {
                        var i = this._segments[e],
                            r = i._handleIn;
                        (i._handleIn = i._handleOut),
                            (i._handleOut = r),
                            (i._index = e);
                    }
                    (this._curves = null),
                        this._clockwise !== t &&
                            (this._clockwise = !this._clockwise),
                        this._changed(9);
                },
                join: function (t) {
                    if (t) {
                        var e = t._segments,
                            n = this.getLastSegment(),
                            i = t.getLastSegment();
                        n._point.equals(i._point) && t.reverse();
                        var r,
                            s = t.getFirstSegment();
                        n._point.equals(s._point)
                            ? (n.setHandleOut(s._handleOut),
                              this._add(e.slice(1)))
                            : ((r = this.getFirstSegment()),
                              r._point.equals(s._point) && t.reverse(),
                              (i = t.getLastSegment()),
                              r._point.equals(i._point)
                                  ? (r.setHandleIn(i._handleIn),
                                    this._add(e.slice(0, e.length - 1), 0))
                                  : this._add(e.slice())),
                            t.closed && this._add([e[0]]),
                            t.remove();
                    }
                    var a = this.getFirstSegment(),
                        o = this.getLastSegment();
                    a !== o &&
                        a._point.equals(o._point) &&
                        (a.setHandleIn(o._handleIn),
                        o.remove(),
                        this.setClosed(!0));
                },
                toShape: function (e) {
                    function n(t, e) {
                        return d[t].isColinear(d[e]);
                    }
                    function i(t) {
                        return d[t].isOrthogonal();
                    }
                    function r(t) {
                        return d[t].isArc();
                    }
                    function s(t, e) {
                        return d[t]._point.getDistance(d[e]._point);
                    }
                    if (!this._closed) return null;
                    var a,
                        h,
                        u,
                        l,
                        d = this._segments;
                    if (
                        (this.isPolygon() &&
                        4 === d.length &&
                        n(0, 2) &&
                        n(1, 3) &&
                        i(1)
                            ? ((a = b.Rectangle),
                              (h = new c(s(0, 3), s(0, 1))),
                              (l = d[1]._point.add(d[2]._point).divide(2)))
                            : 8 === d.length &&
                              r(0) &&
                              r(2) &&
                              r(4) &&
                              r(6) &&
                              n(1, 5) &&
                              n(3, 7)
                            ? ((a = b.Rectangle),
                              (h = new c(s(1, 6), s(0, 3))),
                              (u = h
                                  .subtract(new c(s(0, 7), s(1, 2)))
                                  .divide(2)),
                              (l = d[3]._point.add(d[4]._point).divide(2)))
                            : 4 === d.length &&
                              r(0) &&
                              r(1) &&
                              r(2) &&
                              r(3) &&
                              (o.isZero(s(0, 2) - s(1, 3))
                                  ? ((a = b.Circle), (u = s(0, 2) / 2))
                                  : ((a = b.Ellipse),
                                    (u = new c(s(2, 0) / 2, s(3, 1) / 2))),
                              (l = d[1]._point)),
                        a)
                    ) {
                        var f = this.getPosition(!0),
                            _ = new a({
                                center: f,
                                size: h,
                                radius: u,
                                insert: !1
                            });
                        return (
                            _.rotate(l.subtract(f).getAngle() + 90),
                            _.setStyle(this._style),
                            (e || e === t) && _.insertAbove(this),
                            _
                        );
                    }
                    return null;
                },
                _hitTestSelf: function (t, e) {
                    function n(e, n) {
                        return t.subtract(e).divide(n).length <= 1;
                    }
                    function i(t, i, r) {
                        if (!e.selected || i.isSelected()) {
                            var s = t._point;
                            if ((i !== s && (i = i.add(s)), n(i, w)))
                                return new P(r, _, {
                                    segment: t,
                                    point: i
                                });
                        }
                    }
                    function r(t, n) {
                        return (
                            ((n || e.segments) && i(t, t._point, "segment")) ||
                            (!n &&
                                e.handles &&
                                (i(t, t._handleIn, "handle-in") ||
                                    i(t, t._handleOut, "handle-out")))
                        );
                    }
                    function s(t) {
                        c.add(t);
                    }
                    function a(e) {
                        if (
                            ("round" !== o || "round" !== u) &&
                            ((c = new O({
                                internal: !0,
                                closed: !0
                            })),
                            m || (e._index > 0 && e._index < v - 1)
                                ? "round" !== o &&
                                  (e._handleIn.isZero() ||
                                      e._handleOut.isZero()) &&
                                  O._addBevelJoin(e, o, S, l, s, !0)
                                : "round" !== u &&
                                  O._addSquareCap(e, u, S, s, !0),
                            !c.isEmpty())
                        ) {
                            var i;
                            return (
                                c.contains(t) ||
                                ((i = c.getNearestLocation(t)) &&
                                    n(i.getPoint(), y))
                            );
                        }
                        return n(e._point, w);
                    }
                    var o,
                        u,
                        l,
                        c,
                        d,
                        f,
                        _ = this,
                        g = this.getStyle(),
                        p = this._segments,
                        v = p.length,
                        m = this._closed,
                        y = e._tolerancePadding,
                        w = y,
                        x = e.stroke && g.hasStroke(),
                        b = e.fill && g.hasFill(),
                        C = e.curves,
                        S = x
                            ? g.getStrokeWidth() / 2
                            : (b && e.tolerance > 0) || C
                            ? 0
                            : null;
                    if (
                        (null !== S &&
                            (S > 0
                                ? ((o = g.getStrokeJoin()),
                                  (u = g.getStrokeCap()),
                                  (l = S * g.getMiterLimit()),
                                  (w = y.add(new h(S, S))))
                                : (o = u = "round")),
                        !e.ends || e.segments || m)
                    ) {
                        if (e.segments || e.handles)
                            for (var k = 0; v > k; k++)
                                if ((f = r(p[k]))) return f;
                    } else if ((f = r(p[0], !0) || r(p[v - 1], !0))) return f;
                    if (null !== S) {
                        if ((d = this.getNearestLocation(t))) {
                            var M = d.getParameter();
                            0 === M || (1 === M && v > 1)
                                ? a(d.getSegment()) || (d = null)
                                : n(d.getPoint(), w) || (d = null);
                        }
                        if (!d && "miter" === o && v > 1)
                            for (var k = 0; v > k; k++) {
                                var z = p[k];
                                if (t.getDistance(z._point) <= l && a(z)) {
                                    d = z.getLocation();
                                    break;
                                }
                            }
                    }
                    return (!d && b && this._contains(t)) || (d && !x && !C)
                        ? new P("fill", this)
                        : d
                        ? new P(x ? "stroke" : "curve", this, {
                              location: d,
                              point: d.getPoint()
                          })
                        : null;
                }
            },
            {
                beans: !1,
                _getOffset: function (t) {
                    var e = t && t.getIndex();
                    if (null != e) {
                        for (var n = this.getCurves(), i = 0, r = 0; e > r; r++)
                            i += n[r].getLength();
                        var s = n[e],
                            a = t.getParameter();
                        return a > 0 && (i += s.getPartLength(0, a)), i;
                    }
                    return null;
                },
                getLocationOf: function () {
                    for (
                        var t = h.read(arguments),
                            e = this.getCurves(),
                            n = 0,
                            i = e.length;
                        i > n;
                        n++
                    ) {
                        var r = e[n].getLocationOf(t);
                        if (r) return r;
                    }
                    return null;
                },
                getOffsetOf: function () {
                    var t = this.getLocationOf.apply(this, arguments);
                    return t ? t.getOffset() : null;
                },
                getLocationAt: function (t, e) {
                    var n = this.getCurves(),
                        i = 0;
                    if (e) {
                        var r = ~~t;
                        return n[r].getLocationAt(t - r, !0);
                    }
                    for (var s = 0, a = n.length; a > s; s++) {
                        var o = i,
                            h = n[s];
                        if (((i += h.getLength()), i > t))
                            return h.getLocationAt(t - o);
                    }
                    return t <= this.getLength()
                        ? new A(n[n.length - 1], 1)
                        : null;
                },
                getPointAt: function (t, e) {
                    var n = this.getLocationAt(t, e);
                    return n && n.getPoint();
                },
                getTangentAt: function (t, e) {
                    var n = this.getLocationAt(t, e);
                    return n && n.getTangent();
                },
                getNormalAt: function (t, e) {
                    var n = this.getLocationAt(t, e);
                    return n && n.getNormal();
                },
                getNearestLocation: function () {
                    for (
                        var t = h.read(arguments),
                            e = this.getCurves(),
                            n = 1 / 0,
                            i = null,
                            r = 0,
                            s = e.length;
                        s > r;
                        r++
                    ) {
                        var a = e[r].getNearestLocation(t);
                        a._distance < n && ((n = a._distance), (i = a));
                    }
                    return i;
                },
                getNearestPoint: function () {
                    return this.getNearestLocation
                        .apply(this, arguments)
                        .getPoint();
                }
            },
            new (function () {
                function t(t, e, n, i) {
                    function r(e) {
                        var n = a[e],
                            i = a[e + 1];
                        (c != n || d != i) &&
                            (t.beginPath(),
                            t.moveTo(c, d),
                            t.lineTo(n, i),
                            t.stroke(),
                            t.beginPath(),
                            t.arc(n, i, s, 0, 2 * Math.PI, !0),
                            t.fill());
                    }
                    for (
                        var s = i / 2, a = Array(6), o = 0, h = e.length;
                        h > o;
                        o++
                    ) {
                        var u = e[o];
                        u._transformCoordinates(n, a, !1);
                        var l = u._selectionState,
                            c = a[0],
                            d = a[1];
                        if (
                            (1 & l && r(2),
                            2 & l && r(4),
                            t.fillRect(c - s, d - s, i, i),
                            !(4 & l))
                        ) {
                            var f = t.fillStyle;
                            (t.fillStyle = "#ffffff"),
                                t.fillRect(c - s + 1, d - s + 1, i - 2, i - 2),
                                (t.fillStyle = f);
                        }
                    }
                }
                function e(t, e, n) {
                    function i(e) {
                        if (n)
                            e._transformCoordinates(n, _, !1),
                                (r = _[0]),
                                (s = _[1]);
                        else {
                            var i = e._point;
                            (r = i._x), (s = i._y);
                        }
                        if (g) t.moveTo(r, s), (g = !1);
                        else {
                            if (n) (h = _[2]), (u = _[3]);
                            else {
                                var d = e._handleIn;
                                (h = r + d._x), (u = s + d._y);
                            }
                            h === r && u === s && l === a && c === o
                                ? t.lineTo(r, s)
                                : t.bezierCurveTo(l, c, h, u, r, s);
                        }
                        if (((a = r), (o = s), n)) (l = _[4]), (c = _[5]);
                        else {
                            var d = e._handleOut;
                            (l = a + d._x), (c = o + d._y);
                        }
                    }
                    for (
                        var r,
                            s,
                            a,
                            o,
                            h,
                            u,
                            l,
                            c,
                            d = e._segments,
                            f = d.length,
                            _ = Array(6),
                            g = !0,
                            p = 0;
                        f > p;
                        p++
                    )
                        i(d[p]);
                    e._closed && f > 0 && i(d[0]);
                }
                return {
                    _draw: function (t, n, i) {
                        function r(t) {
                            return l[((t % c) + c) % c];
                        }
                        var s = n.dontStart,
                            a = n.dontFinish || n.clip,
                            o = this.getStyle(),
                            h = o.hasFill(),
                            u = o.hasStroke(),
                            l = o.getDashArray(),
                            c = !paper.support.nativeDash && u && l && l.length;
                        if (
                            (s || t.beginPath(),
                            !s && this._currentPath
                                ? (t.currentPath = this._currentPath)
                                : (h || (u && !c) || a) &&
                                  (e(t, this, i),
                                  this._closed && t.closePath(),
                                  s || (this._currentPath = t.currentPath)),
                            !a &&
                                (h || u) &&
                                (this._setStyles(t),
                                h &&
                                    (t.fill(o.getWindingRule()),
                                    (t.shadowColor = "rgba(0,0,0,0)")),
                                u))
                        ) {
                            if (c) {
                                s || t.beginPath();
                                var d,
                                    f = new L(this, i),
                                    _ = f.length,
                                    g = -o.getDashOffset(),
                                    p = 0;
                                for (g %= _; g > 0; ) g -= r(p--) + r(p--);
                                for (; _ > g; )
                                    (d = g + r(p++)),
                                        (g > 0 || d > 0) &&
                                            f.drawPart(
                                                t,
                                                Math.max(g, 0),
                                                Math.max(d, 0)
                                            ),
                                        (g = d + r(p++));
                            }
                            t.stroke();
                        }
                    },
                    _drawSelected: function (n, i) {
                        n.beginPath(),
                            e(n, this, i),
                            n.stroke(),
                            t(n, this._segments, i, paper.settings.handleSize);
                    }
                };
            })(),
            new (function () {
                function t(t) {
                    var e = t.length,
                        n = [],
                        i = [],
                        r = 2;
                    n[0] = t[0] / r;
                    for (var s = 1; e > s; s++)
                        (i[s] = 1 / r),
                            (r = (e - 1 > s ? 4 : 2) - i[s]),
                            (n[s] = (t[s] - n[s - 1]) / r);
                    for (var s = 1; e > s; s++)
                        n[e - s - 1] -= i[e - s] * n[e - s];
                    return n;
                }
                return {
                    smooth: function () {
                        var e = this._segments,
                            n = e.length,
                            i = this._closed,
                            r = n,
                            s = 0;
                        if (!(2 >= n)) {
                            i &&
                                ((s = Math.min(n, 4)),
                                (r += 2 * Math.min(n, s)));
                            for (var a = [], o = 0; n > o; o++)
                                a[o + s] = e[o]._point;
                            if (i)
                                for (var o = 0; s > o; o++)
                                    (a[o] = e[o + n - s]._point),
                                        (a[o + n + s] = e[o]._point);
                            else r--;
                            for (var u = [], o = 1; r - 1 > o; o++)
                                u[o] = 4 * a[o]._x + 2 * a[o + 1]._x;
                            (u[0] = a[0]._x + 2 * a[1]._x),
                                (u[r - 1] = 3 * a[r - 1]._x);
                            for (var l = t(u), o = 1; r - 1 > o; o++)
                                u[o] = 4 * a[o]._y + 2 * a[o + 1]._y;
                            (u[0] = a[0]._y + 2 * a[1]._y),
                                (u[r - 1] = 3 * a[r - 1]._y);
                            var c = t(u);
                            if (i) {
                                for (var o = 0, d = n; s > o; o++, d++) {
                                    var f = o / s,
                                        _ = 1 - f,
                                        g = o + s,
                                        p = d + s;
                                    (l[d] = l[o] * f + l[d] * _),
                                        (c[d] = c[o] * f + c[d] * _),
                                        (l[p] = l[g] * _ + l[p] * f),
                                        (c[p] = c[g] * _ + c[p] * f);
                                }
                                r--;
                            }
                            for (var v = null, o = s; r - s >= o; o++) {
                                var m = e[o - s];
                                v && m.setHandleIn(v.subtract(m._point)),
                                    r > o &&
                                        (m.setHandleOut(
                                            new h(l[o], c[o]).subtract(m._point)
                                        ),
                                        (v =
                                            r - 1 > o
                                                ? new h(
                                                      2 * a[o + 1]._x -
                                                          l[o + 1],
                                                      2 * a[o + 1]._y - c[o + 1]
                                                  )
                                                : new h(
                                                      (a[r]._x + l[r - 1]) / 2,
                                                      (a[r]._y + c[r - 1]) / 2
                                                  )));
                            }
                            if (i && v) {
                                var m = this._segments[0];
                                m.setHandleIn(v.subtract(m._point));
                            }
                        }
                    }
                };
            })(),
            new (function () {
                function t(t) {
                    var e = t._segments;
                    if (0 === e.length)
                        throw Error("Use a moveTo() command first");
                    return e[e.length - 1];
                }
                return {
                    moveTo: function () {
                        var t = this._segments;
                        1 === t.length && this.removeSegment(0),
                            t.length || this._add([new k(h.read(arguments))]);
                    },
                    moveBy: function () {
                        throw Error("moveBy() is unsupported on Path items.");
                    },
                    lineTo: function () {
                        this._add([new k(h.read(arguments))]);
                    },
                    cubicCurveTo: function () {
                        var e = h.read(arguments),
                            n = h.read(arguments),
                            i = h.read(arguments),
                            r = t(this);
                        r.setHandleOut(e.subtract(r._point)),
                            this._add([new k(i, n.subtract(i))]);
                    },
                    quadraticCurveTo: function () {
                        var e = h.read(arguments),
                            n = h.read(arguments),
                            i = t(this)._point;
                        this.cubicCurveTo(
                            e.add(i.subtract(e).multiply(1 / 3)),
                            e.add(n.subtract(e).multiply(1 / 3)),
                            n
                        );
                    },
                    curveTo: function () {
                        var n = h.read(arguments),
                            i = h.read(arguments),
                            r = e.pick(e.read(arguments), 0.5),
                            s = 1 - r,
                            a = t(this)._point,
                            o = n
                                .subtract(a.multiply(s * s))
                                .subtract(i.multiply(r * r))
                                .divide(2 * r * s);
                        if (o.isNaN())
                            throw Error(
                                "Cannot put a curve through points with parameter = " +
                                    r
                            );
                        this.quadraticCurveTo(o, i);
                    },
                    arcTo: function () {
                        var n,
                            i,
                            r,
                            s,
                            a,
                            o = t(this),
                            u = o._point,
                            l = h.read(arguments),
                            d = e.peek(arguments),
                            f = e.pick(d, !0);
                        if ("boolean" == typeof f)
                            var _ = u.add(l).divide(2),
                                n = _.add(_.subtract(u).rotate(f ? -90 : 90));
                        else if (e.remain(arguments) <= 2)
                            (n = l), (l = h.read(arguments));
                        else {
                            var v = c.read(arguments);
                            if (v.isZero()) return this.lineTo(l);
                            var m = e.read(arguments),
                                f = !!e.read(arguments),
                                y = !!e.read(arguments),
                                _ = u.add(l).divide(2),
                                w = u.subtract(_).rotate(-m),
                                x = w.x,
                                b = w.y,
                                C = Math.abs,
                                S = 1e-11,
                                P = C(v.width),
                                M = C(v.height),
                                z = P * P,
                                A = M * M,
                                I = x * x,
                                O = b * b,
                                T = Math.sqrt(I / z + O / A);
                            if (
                                (T > 1 &&
                                    ((P *= T),
                                    (M *= T),
                                    (z = P * P),
                                    (A = M * M)),
                                (T = (z * A - z * O - A * I) / (z * O + A * I)),
                                C(T) < S && (T = 0),
                                0 > T)
                            )
                                throw Error(
                                    "Cannot create an arc with the given arguments"
                                );
                            (i = new h((P * b) / M, (-M * x) / P)
                                .multiply((y === f ? -1 : 1) * Math.sqrt(T))
                                .rotate(m)
                                .add(_)),
                                (a = new g()
                                    .translate(i)
                                    .rotate(m)
                                    .scale(P, M)),
                                (s = a._inverseTransform(u)),
                                (r = s.getDirectedAngle(
                                    a._inverseTransform(l)
                                )),
                                !f && r > 0
                                    ? (r -= 360)
                                    : f && 0 > r && (r += 360);
                        }
                        if (n) {
                            var L = new p(
                                    u.add(n).divide(2),
                                    n.subtract(u).rotate(90),
                                    !0
                                ),
                                E = new p(
                                    n.add(l).divide(2),
                                    l.subtract(n).rotate(90),
                                    !0
                                ),
                                N = new p(u, l),
                                j = N.getSide(n);
                            if (((i = L.intersect(E, !0)), !i)) {
                                if (!j) return this.lineTo(l);
                                throw Error(
                                    "Cannot create an arc with the given arguments"
                                );
                            }
                            (s = u.subtract(i)),
                                (r = s.getDirectedAngle(l.subtract(i)));
                            var D = N.getSide(i);
                            0 === D
                                ? (r = j * Math.abs(r))
                                : j === D && (r += 0 > r ? 360 : -360);
                        }
                        for (
                            var B = Math.abs(r),
                                R = B >= 360 ? 4 : Math.ceil(B / 90),
                                F = r / R,
                                V = (F * Math.PI) / 360,
                                q = ((4 / 3) * Math.sin(V)) / (1 + Math.cos(V)),
                                Z = [],
                                U = 0;
                            R >= U;
                            U++
                        ) {
                            var w = l,
                                H = null;
                            if (
                                (R > U &&
                                    ((H = s.rotate(90).multiply(q)),
                                    a
                                        ? ((w = a._transformPoint(s)),
                                          (H = a
                                              ._transformPoint(s.add(H))
                                              .subtract(w)))
                                        : (w = i.add(s))),
                                0 === U)
                            )
                                o.setHandleOut(H);
                            else {
                                var W = s.rotate(-90).multiply(q);
                                a &&
                                    (W = a
                                        ._transformPoint(s.add(W))
                                        .subtract(w)),
                                    Z.push(new k(w, W, H));
                            }
                            s = s.rotate(F);
                        }
                        this._add(Z);
                    },
                    lineBy: function () {
                        var e = h.read(arguments),
                            n = t(this)._point;
                        this.lineTo(n.add(e));
                    },
                    curveBy: function () {
                        var n = h.read(arguments),
                            i = h.read(arguments),
                            r = e.read(arguments),
                            s = t(this)._point;
                        this.curveTo(s.add(n), s.add(i), r);
                    },
                    cubicCurveBy: function () {
                        var e = h.read(arguments),
                            n = h.read(arguments),
                            i = h.read(arguments),
                            r = t(this)._point;
                        this.cubicCurveTo(r.add(e), r.add(n), r.add(i));
                    },
                    quadraticCurveBy: function () {
                        var e = h.read(arguments),
                            n = h.read(arguments),
                            i = t(this)._point;
                        this.quadraticCurveTo(i.add(e), i.add(n));
                    },
                    arcBy: function () {
                        var n = t(this)._point,
                            i = n.add(h.read(arguments)),
                            r = e.pick(e.peek(arguments), !0);
                        "boolean" == typeof r
                            ? this.arcTo(i, r)
                            : this.arcTo(i, n.add(h.read(arguments)));
                    },
                    closePath: function (t) {
                        this.setClosed(!0), t && this.join();
                    }
                };
            })(),
            {
                _getBounds: function (t, e) {
                    return O[t](
                        this._segments,
                        this._closed,
                        this.getStyle(),
                        e
                    );
                },
                statics: {
                    isClockwise: function (t) {
                        for (var e = 0, n = 0, i = t.length; i > n; n++)
                            for (
                                var r = z.getValues(
                                        t[n],
                                        t[i > n + 1 ? n + 1 : 0]
                                    ),
                                    s = 2;
                                8 > s;
                                s += 2
                            )
                                e += (r[s - 2] - r[s]) * (r[s + 1] + r[s - 1]);
                        return e > 0;
                    },
                    getBounds: function (t, e, n, i, r) {
                        function s(t) {
                            t._transformCoordinates(i, o, !1);
                            for (var e = 0; 2 > e; e++)
                                z._addBounds(
                                    h[e],
                                    h[e + 4],
                                    o[e + 2],
                                    o[e],
                                    e,
                                    r ? r[e] : 0,
                                    u,
                                    l,
                                    c
                                );
                            var n = h;
                            (h = o), (o = n);
                        }
                        var a = t[0];
                        if (!a) return new f();
                        for (
                            var o = Array(6),
                                h = a._transformCoordinates(i, Array(6), !1),
                                u = h.slice(0, 2),
                                l = u.slice(),
                                c = Array(2),
                                d = 1,
                                _ = t.length;
                            _ > d;
                            d++
                        )
                            s(t[d]);
                        return (
                            e && s(a),
                            new f(u[0], u[1], l[0] - u[0], l[1] - u[1])
                        );
                    },
                    getStrokeBounds: function (t, e, n, i) {
                        function r(t) {
                            d = d.include(i ? i._transformPoint(t, t) : t);
                        }
                        function s(t) {
                            d = d.unite(
                                v.setCenter(
                                    i ? i._transformPoint(t._point) : t._point
                                )
                            );
                        }
                        function a(t, e) {
                            var n = t._handleIn,
                                i = t._handleOut;
                            "round" === e ||
                            (!n.isZero() && !i.isZero() && n.isColinear(i))
                                ? s(t)
                                : O._addBevelJoin(t, e, u, p, r);
                        }
                        function o(t, e) {
                            "round" === e ? s(t) : O._addSquareCap(t, e, u, r);
                        }
                        if (!n.hasStroke()) return O.getBounds(t, e, n, i);
                        for (
                            var h = t.length - (e ? 0 : 1),
                                u = n.getStrokeWidth() / 2,
                                l = O._getPenPadding(u, i),
                                d = O.getBounds(t, e, n, i, l),
                                _ = n.getStrokeJoin(),
                                g = n.getStrokeCap(),
                                p = u * n.getMiterLimit(),
                                v = new f(new c(l).multiply(2)),
                                m = 1;
                            h > m;
                            m++
                        )
                            a(t[m], _);
                        return (
                            e
                                ? a(t[0], _)
                                : h > 0 && (o(t[0], g), o(t[t.length - 1], g)),
                            d
                        );
                    },
                    _getPenPadding: function (t, e) {
                        if (!e) return [t, t];
                        var n = e.shiftless(),
                            i = n.transform(new h(t, 0)),
                            r = n.transform(new h(0, t)),
                            s = i.getAngleInRadians(),
                            a = i.getLength(),
                            o = r.getLength(),
                            u = Math.sin(s),
                            l = Math.cos(s),
                            c = Math.tan(s),
                            d = -Math.atan((o * c) / a),
                            f = Math.atan(o / (c * a));
                        return [
                            Math.abs(a * Math.cos(d) * l - o * Math.sin(d) * u),
                            Math.abs(o * Math.sin(f) * l + a * Math.cos(f) * u)
                        ];
                    },
                    _addBevelJoin: function (t, e, n, i, r, s) {
                        var a = t.getCurve(),
                            o = a.getPrevious(),
                            u = a.getPointAt(0, !0),
                            l = o.getNormalAt(1, !0),
                            c = a.getNormalAt(0, !0),
                            d = l.getDirectedAngle(c) < 0 ? -n : n;
                        if (
                            (l.setLength(d),
                            c.setLength(d),
                            s && (r(u), r(u.add(l))),
                            "miter" === e)
                        ) {
                            var f = new p(
                                u.add(l),
                                new h(-l.y, l.x),
                                !0
                            ).intersect(
                                new p(u.add(c), new h(-c.y, c.x), !0),
                                !0
                            );
                            if (f && u.getDistance(f) <= i && (r(f), !s))
                                return;
                        }
                        s || r(u.add(l)), r(u.add(c));
                    },
                    _addSquareCap: function (t, e, n, i, r) {
                        var s = t._point,
                            a = t.getLocation(),
                            o = a.getNormal().normalize(n);
                        r && (i(s.subtract(o)), i(s.add(o))),
                            "square" === e &&
                                (s = s.add(
                                    o.rotate(0 === a.getParameter() ? -90 : 90)
                                )),
                            i(s.add(o)),
                            i(s.subtract(o));
                    },
                    getHandleBounds: function (t, e, n, i, r, s) {
                        for (
                            var a = Array(6),
                                o = 1 / 0,
                                h = -o,
                                u = o,
                                l = h,
                                c = 0,
                                d = t.length;
                            d > c;
                            c++
                        ) {
                            var _ = t[c];
                            _._transformCoordinates(i, a, !1);
                            for (var g = 0; 6 > g; g += 2) {
                                var p = 0 === g ? s : r,
                                    v = p ? p[0] : 0,
                                    m = p ? p[1] : 0,
                                    y = a[g],
                                    w = a[g + 1],
                                    x = y - v,
                                    b = y + v,
                                    C = w - m,
                                    S = w + m;
                                o > x && (o = x),
                                    b > h && (h = b),
                                    u > C && (u = C),
                                    S > l && (l = S);
                            }
                        }
                        return new f(o, u, h - o, l - u);
                    },
                    getRoughBounds: function (t, e, n, i) {
                        var r = n.hasStroke() ? n.getStrokeWidth() / 2 : 0,
                            s = r;
                        return (
                            r > 0 &&
                                ("miter" === n.getStrokeJoin() &&
                                    (s = r * n.getMiterLimit()),
                                "square" === n.getStrokeCap() &&
                                    (s = Math.max(s, r * Math.sqrt(2)))),
                            O.getHandleBounds(
                                t,
                                e,
                                n,
                                i,
                                O._getPenPadding(r, i),
                                O._getPenPadding(s, i)
                            )
                        );
                    }
                }
            }
        );
    O.inject({
        statics: new (function () {
            function t(t, n, i) {
                var r = e.getNamed(i),
                    s = new O(r && r.insert === !1 && y.NO_INSERT);
                return s._add(t), (s._closed = n), s.set(r);
            }
            function n(e, n, i) {
                for (var s = Array(4), a = 0; 4 > a; a++) {
                    var o = r[a];
                    s[a] = new k(
                        o._point.multiply(n).add(e),
                        o._handleIn.multiply(n),
                        o._handleOut.multiply(n)
                    );
                }
                return t(s, !0, i);
            }
            var i = 0.5522847498307936,
                r = [
                    new k([-1, 0], [0, i], [0, -i]),
                    new k([0, -1], [-i, 0], [i, 0]),
                    new k([1, 0], [0, -i], [0, i]),
                    new k([0, 1], [i, 0], [-i, 0])
                ];
            return {
                Line: function () {
                    return t(
                        [
                            new k(h.readNamed(arguments, "from")),
                            new k(h.readNamed(arguments, "to"))
                        ],
                        !1,
                        arguments
                    );
                },
                Circle: function () {
                    var t = h.readNamed(arguments, "center"),
                        i = e.readNamed(arguments, "radius");
                    return n(t, new c(i), arguments);
                },
                Rectangle: function () {
                    var e,
                        n = f.readNamed(arguments, "rectangle"),
                        r = c.readNamed(arguments, "radius", 0, {
                            readNull: !0
                        }),
                        s = n.getBottomLeft(!0),
                        a = n.getTopLeft(!0),
                        o = n.getTopRight(!0),
                        h = n.getBottomRight(!0);
                    if (!r || r.isZero())
                        e = [new k(s), new k(a), new k(o), new k(h)];
                    else {
                        r = c.min(r, n.getSize(!0).divide(2));
                        var u = r.width,
                            l = r.height,
                            d = u * i,
                            _ = l * i;
                        e = [
                            new k(s.add(u, 0), null, [-d, 0]),
                            new k(s.subtract(0, l), [0, _]),
                            new k(a.add(0, l), null, [0, -_]),
                            new k(a.add(u, 0), [-d, 0], null),
                            new k(o.subtract(u, 0), null, [d, 0]),
                            new k(o.add(0, l), [0, -_], null),
                            new k(h.subtract(0, l), null, [0, _]),
                            new k(h.subtract(u, 0), [d, 0])
                        ];
                    }
                    return t(e, !0, arguments);
                },
                RoundRectangle: "#Rectangle",
                Ellipse: function () {
                    var t = b._readEllipse(arguments);
                    return n(t.center, t.radius, arguments);
                },
                Oval: "#Ellipse",
                Arc: function () {
                    var t = h.readNamed(arguments, "from"),
                        n = h.readNamed(arguments, "through"),
                        i = h.readNamed(arguments, "to"),
                        r = e.getNamed(arguments),
                        s = new O(r && r.insert === !1 && y.NO_INSERT);
                    return s.moveTo(t), s.arcTo(n, i), s.set(r);
                },
                RegularPolygon: function () {
                    for (
                        var n = h.readNamed(arguments, "center"),
                            i = e.readNamed(arguments, "sides"),
                            r = e.readNamed(arguments, "radius"),
                            s = 360 / i,
                            a = !(i % 3),
                            o = new h(0, a ? -r : r),
                            u = a ? -1 : 0.5,
                            l = Array(i),
                            c = 0;
                        i > c;
                        c++
                    )
                        l[c] = new k(n.add(o.rotate((c + u) * s)));
                    return t(l, !0, arguments);
                },
                Star: function () {
                    for (
                        var n = h.readNamed(arguments, "center"),
                            i = 2 * e.readNamed(arguments, "points"),
                            r = e.readNamed(arguments, "radius1"),
                            s = e.readNamed(arguments, "radius2"),
                            a = 360 / i,
                            o = new h(0, -1),
                            u = Array(i),
                            l = 0;
                        i > l;
                        l++
                    )
                        u[l] = new k(
                            n.add(o.rotate(a * l).multiply(l % 2 ? s : r))
                        );
                    return t(u, !0, arguments);
                }
            };
        })()
    });
    var T = I.extend(
        {
            _class: "CompoundPath",
            _serializeFields: {
                children: []
            },
            initialize: function (t) {
                (this._children = []),
                    (this._namedChildren = {}),
                    this._initialize(t) ||
                        ("string" == typeof t
                            ? this.setPathData(t)
                            : this.addChildren(
                                  Array.isArray(t) ? t : arguments
                              ));
            },
            insertChildren: function xe(e, n, i) {
                n = xe.base.call(this, e, n, i, O);
                for (var r = 0, s = !i && n && n.length; s > r; r++) {
                    var a = n[r];
                    a._clockwise === t && a.setClockwise(0 === a._index);
                }
                return n;
            },
            reverse: function () {
                for (var t = this._children, e = 0, n = t.length; n > e; e++)
                    t[e].reverse();
            },
            smooth: function () {
                for (var t = 0, e = this._children.length; e > t; t++)
                    this._children[t].smooth();
            },
            isClockwise: function () {
                var t = this.getFirstChild();
                return t && t.isClockwise();
            },
            setClockwise: function (t) {
                this.isClockwise() !== !!t && this.reverse();
            },
            getFirstSegment: function () {
                var t = this.getFirstChild();
                return t && t.getFirstSegment();
            },
            getLastSegment: function () {
                var t = this.getLastChild();
                return t && t.getLastSegment();
            },
            getCurves: function () {
                for (
                    var t = this._children, e = [], n = 0, i = t.length;
                    i > n;
                    n++
                )
                    e.push.apply(e, t[n].getCurves());
                return e;
            },
            getFirstCurve: function () {
                var t = this.getFirstChild();
                return t && t.getFirstCurve();
            },
            getLastCurve: function () {
                var t = this.getLastChild();
                return t && t.getFirstCurve();
            },
            getArea: function () {
                for (
                    var t = this._children, e = 0, n = 0, i = t.length;
                    i > n;
                    n++
                )
                    e += t[n].getArea();
                return e;
            }
        },
        {
            beans: !0,
            getPathData: function (t, e) {
                for (
                    var n = this._children, i = [], r = 0, s = n.length;
                    s > r;
                    r++
                ) {
                    var a = n[r],
                        o = a._matrix;
                    i.push(
                        a.getPathData(t && !o.isIdentity() ? t.chain(o) : o, e)
                    );
                }
                return i.join(" ");
            }
        },
        {
            _getChildHitTestOptions: function (t) {
                return t.class === O || "path" === t.type
                    ? t
                    : new e(t, {
                          fill: !1
                      });
            },
            _draw: function (t, e, n) {
                var i = this._children;
                if (0 !== i.length) {
                    if (this._currentPath) t.currentPath = this._currentPath;
                    else {
                        (e = e.extend({
                            dontStart: !0,
                            dontFinish: !0
                        })),
                            t.beginPath();
                        for (var r = 0, s = i.length; s > r; r++)
                            i[r].draw(t, e, n);
                        this._currentPath = t.currentPath;
                    }
                    if (!e.clip) {
                        this._setStyles(t);
                        var a = this._style;
                        a.hasFill() &&
                            (t.fill(a.getWindingRule()),
                            (t.shadowColor = "rgba(0,0,0,0)")),
                            a.hasStroke() && t.stroke();
                    }
                }
            },
            _drawSelected: function (t, e, n) {
                for (var i = this._children, r = 0, s = i.length; s > r; r++) {
                    var a = i[r],
                        o = a._matrix;
                    n[a._id] ||
                        a._drawSelected(t, o.isIdentity() ? e : e.chain(o));
                }
            }
        },
        new (function () {
            function t(t, e) {
                var n = t._children;
                if (e && 0 === n.length)
                    throw Error("Use a moveTo() command first");
                return n[n.length - 1];
            }
            var n = {
                moveTo: function () {
                    var e = t(this),
                        n = e && e.isEmpty() ? e : new O();
                    n !== e && this.addChild(n), n.moveTo.apply(n, arguments);
                },
                moveBy: function () {
                    var e = t(this, !0),
                        n = e && e.getLastSegment(),
                        i = h.read(arguments);
                    this.moveTo(n ? i.add(n._point) : i);
                },
                closePath: function (e) {
                    t(this, !0).closePath(e);
                }
            };
            return (
                e.each(
                    [
                        "lineTo",
                        "cubicCurveTo",
                        "quadraticCurveTo",
                        "curveTo",
                        "arcTo",
                        "lineBy",
                        "cubicCurveBy",
                        "quadraticCurveBy",
                        "curveBy",
                        "arcBy"
                    ],
                    function (e) {
                        n[e] = function () {
                            var n = t(this, !0);
                            n[e].apply(n, arguments);
                        };
                    }
                ),
                n
            );
        })()
    );
    I.inject(
        new (function () {
            function t(t, r, s, a) {
                function o(t) {
                    return t.clone(!1).reduce().reorient().transform(null, !0);
                }
                function h(t) {
                    for (var e = 0, n = t.length; n > e; e++) {
                        var i = t[e];
                        _.push.apply(_, i._segments),
                            g.push.apply(g, i._getMonoCurves());
                    }
                }
                var u = o(t),
                    l = r && t !== r && o(r);
                u.isClockwise() || u.reverse(),
                    !l || a ^ l.isClockwise() || l.reverse(),
                    e(u.getIntersections(l, !0));
                var c = [],
                    d = [],
                    f = [],
                    _ = [],
                    g = [];
                h(u._children || [u]),
                    l && h(l._children || [l]),
                    _.sort(function (t, e) {
                        var n = t._intersection,
                            i = e._intersection;
                        return (!n && !i) || (n && i) ? 0 : n ? -1 : 1;
                    });
                for (var p = 0, v = _.length; v > p; p++) {
                    var m = _[p];
                    if (null == m._winding) {
                        c.length = d.length = f.length = 0;
                        var y = 0,
                            w = m;
                        do
                            c.push(m),
                                f.push((y += m.getCurve().getLength())),
                                (m = m.getNext());
                        while (m && !m._intersection && m !== w);
                        for (var x = 0; 3 > x; x++) {
                            var b = y * Math.random(),
                                C = f.length,
                                S = 0;
                            do
                                if (f[S] >= b) {
                                    S > 0 && (b -= f[S - 1]);
                                    break;
                                }
                            while (++S < C);
                            var P = c[S].getCurve(),
                                k = P.getPointAt(b),
                                M = P.isHorizontal(),
                                z = P._path;
                            z._parent instanceof T && (z = z._parent),
                                (d[x] =
                                    a &&
                                    l &&
                                    ((z === u && l._getWinding(k, M)) ||
                                        (z === l && !u._getWinding(k, M)))
                                        ? 0
                                        : n(k, g, M));
                        }
                        d.sort();
                        for (var A = d[1], x = c.length - 1; x >= 0; x--)
                            c[x]._winding = A;
                    }
                }
                var I = new T();
                return (
                    I.addChildren(i(_, s), !0),
                    u.remove(),
                    l && l.remove(),
                    (I = I.reduce()),
                    I.setStyle(t._style),
                    I
                );
            }
            function e(t) {
                function e() {
                    for (var t = 0, e = n.length; e > t; t++) {
                        var i = n[t];
                        i._handleOut.set(0, 0), i._handleIn.set(0, 0);
                    }
                }
                for (var n, i, r, s = 1e-5, a = t.length - 1; a >= 0; a--) {
                    var o = t[a],
                        h = o._parameter;
                    r && r._curve === o._curve && r._parameter > 0
                        ? (h /= r._parameter)
                        : (n && e(), (i = o._curve), (n = i.isLinear() && []));
                    var u, l;
                    (u = i.divide(h, !0, !0))
                        ? ((l = u._segment1), (i = u.getPrevious()))
                        : (l =
                              s > h
                                  ? i._segment1
                                  : h > 1 - s
                                  ? i._segment2
                                  : i.getPartLength(0, h) <
                                    i.getPartLength(h, 1)
                                  ? i._segment1
                                  : i._segment2),
                        (l._intersection = o.getIntersection()),
                        (o._segment = l),
                        n && n.push(l),
                        (r = o);
                }
                n && e();
            }
            function n(t, e, i, r) {
                var s = 1e-5,
                    a = t.x,
                    o = t.y,
                    u = 0,
                    l = 0,
                    c = [],
                    d = Math.abs,
                    f = 1 - s;
                if (i) {
                    for (
                        var _ = -1 / 0,
                            g = 1 / 0,
                            p = o - s,
                            v = o + s,
                            m = 0,
                            y = e.length;
                        y > m;
                        m++
                    ) {
                        var w = e[m].values;
                        if (z.solveCubic(w, 0, a, c, 0, 1) > 0)
                            for (var x = c.length - 1; x >= 0; x--) {
                                var b = z.evaluate(w, c[x], 0).y;
                                p > b && b > _
                                    ? (_ = b)
                                    : b > v && g > b && (g = b);
                            }
                    }
                    (_ = (_ + o) / 2),
                        (g = (g + o) / 2),
                        _ > -1 / 0 && (u = n(new h(a, _), e)),
                        1 / 0 > g && (l = n(new h(a, g), e));
                } else
                    for (
                        var C = a - s, S = a + s, m = 0, y = e.length;
                        y > m;
                        m++
                    ) {
                        var P = e[m],
                            w = P.values,
                            k = P.winding,
                            M = P.next;
                        if (
                            k &&
                            ((1 === k && o >= w[1] && o <= w[7]) ||
                                (o >= w[7] && o <= w[1])) &&
                            1 ===
                                z.solveCubic(
                                    w,
                                    1,
                                    o,
                                    c,
                                    0,
                                    M.winding || M.values[1] !== o ? f : 1
                                )
                        ) {
                            var A = c[0],
                                I = z.evaluate(w, A, 0).x,
                                O = z.evaluate(w, A, 1).y;
                            (d(O) < s && !z.isLinear(w)) ||
                            (s > A &&
                                O * z.evaluate(P.previous.values, A, 1).y < 0)
                                ? r && I >= C && S >= I && (++u, ++l)
                                : C >= I
                                ? (u += k)
                                : I >= S && (l += k);
                        }
                    }
                return Math.max(d(u), d(l));
            }
            function i(t, e, n) {
                e =
                    e ||
                    function () {
                        return !0;
                    };
                for (
                    var i, r, s = [], a = 0.001, o = 0.999, h = 0, u = t.length;
                    u > h;
                    h++
                )
                    if (((i = r = t[h]), !i._visited && e(i._winding))) {
                        var l = new O(y.NO_INSERT),
                            c = i._intersection,
                            d = c && c._segment,
                            f = !1,
                            _ = 1;
                        do {
                            var g,
                                p = _ > 0 ? i._handleIn : i._handleOut,
                                v = _ > 0 ? i._handleOut : i._handleIn;
                            if (
                                f &&
                                (!e(i._winding) || n) &&
                                (c = i._intersection) &&
                                (g = c._segment) &&
                                g !== r
                            ) {
                                if (n)
                                    (i._visited = g._visited), (i = g), (_ = 1);
                                else {
                                    var m = i.getCurve();
                                    _ > 0 && (m = m.getPrevious());
                                    var w = m.getTangentAt(1 > _ ? a : o, !0),
                                        x = g.getCurve(),
                                        b = x.getPrevious(),
                                        C = b.getTangentAt(o, !0),
                                        S = x.getTangentAt(a, !0),
                                        P = w.cross(C),
                                        M = w.cross(S);
                                    if (0 !== P * M) {
                                        var z = M > P ? b : x,
                                            A = e(z._segment1._winding)
                                                ? z
                                                : M > P
                                                ? x
                                                : b,
                                            I = A._segment1;
                                        (_ = A === b ? -1 : 1),
                                            (I._visited &&
                                                i._path !== I._path) ||
                                            !e(I._winding)
                                                ? (_ = 1)
                                                : ((i._visited = g._visited),
                                                  (i = g),
                                                  I._visited && (_ = 1));
                                    } else _ = 1;
                                }
                                v = _ > 0 ? i._handleOut : i._handleIn;
                            }
                            l.add(new k(i._point, f && p, v)),
                                (f = !0),
                                (i._visited = !0),
                                (i = _ > 0 ? i.getNext() : i.getPrevious());
                        } while (
                            i &&
                            !i._visited &&
                            i !== r &&
                            i !== d &&
                            (i._intersection || e(i._winding))
                        );
                        !i || (i !== r && i !== d)
                            ? l.lastSegment._handleOut.set(0, 0)
                            : (l.firstSegment.setHandleIn(
                                  (i === d ? d : i)._handleIn
                              ),
                              l.setClosed(!0)),
                            l._segments.length >
                                (l._closed ? (l.isPolygon() ? 2 : 0) : 1) &&
                                s.push(l);
                    }
                return s;
            }
            return {
                _getWinding: function (t, e, i) {
                    return n(t, this._getMonoCurves(), e, i);
                },
                unite: function (e) {
                    return t(
                        this,
                        e,
                        function (t) {
                            return 1 === t || 0 === t;
                        },
                        !1
                    );
                },
                intersect: function (e) {
                    return t(
                        this,
                        e,
                        function (t) {
                            return 2 === t;
                        },
                        !1
                    );
                },
                subtract: function (e) {
                    return t(
                        this,
                        e,
                        function (t) {
                            return 1 === t;
                        },
                        !0
                    );
                },
                exclude: function (t) {
                    return new w([this.subtract(t), t.subtract(this)]);
                },
                divide: function (t) {
                    return new w([this.subtract(t), this.intersect(t)]);
                }
            };
        })()
    ),
        O.inject({
            _getMonoCurves: function () {
                function t(t) {
                    var e = t[1],
                        r = t[7],
                        s = {
                            values: t,
                            winding: e === r ? 0 : e > r ? -1 : 1,
                            previous: n,
                            next: null
                        };
                    n && (n.next = s), i.push(s), (n = s);
                }
                function e(e) {
                    if (0 !== z.getLength(e)) {
                        var n = e[1],
                            i = e[3],
                            r = e[5],
                            s = e[7];
                        if (z.isLinear(e)) t(e);
                        else {
                            var a = 3 * (i - r) - n + s,
                                h = 2 * (n + r) - 4 * i,
                                u = i - n,
                                l = 1e-5,
                                c = [],
                                d = o.solveQuadratic(a, h, u, c, l, 1 - l);
                            if (0 === d) t(e);
                            else {
                                c.sort();
                                var f = c[0],
                                    _ = z.subdivide(e, f);
                                t(_[0]),
                                    d > 1 &&
                                        ((f = (c[1] - f) / (1 - f)),
                                        (_ = z.subdivide(_[1], f)),
                                        t(_[0])),
                                    t(_[1]);
                            }
                        }
                    }
                }
                var n,
                    i = this._monoCurves;
                if (!i) {
                    i = this._monoCurves = [];
                    for (
                        var r = this.getCurves(),
                            s = this._segments,
                            a = 0,
                            h = r.length;
                        h > a;
                        a++
                    )
                        e(r[a].getValues());
                    if (!this._closed && s.length > 1) {
                        var u = s[s.length - 1]._point,
                            l = s[0]._point,
                            c = u._x,
                            d = u._y,
                            f = l._x,
                            _ = l._y;
                        e([c, d, c, d, f, _, f, _]);
                    }
                    if (i.length > 0) {
                        var g = i[0],
                            p = i[i.length - 1];
                        (g.previous = p), (p.next = g);
                    }
                }
                return i;
            },
            getInteriorPoint: function () {
                var t = this.getBounds(),
                    e = t.getCenter(!0);
                if (!this.contains(e)) {
                    for (
                        var n = this._getMonoCurves(),
                            i = [],
                            r = e.y,
                            s = [],
                            a = 0,
                            o = n.length;
                        o > a;
                        a++
                    ) {
                        var h = n[a].values;
                        if (
                            ((1 === n[a].winding && r >= h[1] && r <= h[7]) ||
                                (r >= h[7] && r <= h[1])) &&
                            z.solveCubic(h, 1, r, i, 0, 1) > 0
                        )
                            for (var u = i.length - 1; u >= 0; u--)
                                s.push(z.evaluate(h, i[u], 0).x);
                        if (s.length > 1) break;
                    }
                    e.x = (s[0] + s[1]) / 2;
                }
                return e;
            },
            reorient: function () {
                return this.setClockwise(!0), this;
            }
        }),
        T.inject({
            _getMonoCurves: function () {
                for (
                    var t = this._children, e = [], n = 0, i = t.length;
                    i > n;
                    n++
                )
                    e.push.apply(e, t[n]._getMonoCurves());
                return e;
            },
            reorient: function () {
                var t = this.removeChildren().sort(function (t, e) {
                    return e.getBounds().getArea() - t.getBounds().getArea();
                });
                this.addChildren(t);
                for (
                    var e = t[0].isClockwise(), n = 1, i = t.length;
                    i > n;
                    n++
                ) {
                    for (
                        var r = t[n].getInteriorPoint(), s = 0, a = n - 1;
                        a >= 0;
                        a--
                    )
                        t[a].contains(r) && s++;
                    t[n].setClockwise(0 === s % 2 && e);
                }
                return this;
            }
        });
    var L = e.extend({
            initialize: function (t, e) {
                function n(t, n) {
                    var i = z.getValues(t, n, e);
                    a.curves.push(i), a._computeParts(i, t._index, 0, 1);
                }
                (this.curves = []),
                    (this.parts = []),
                    (this.length = 0),
                    (this.index = 0);
                for (
                    var i,
                        r = t._segments,
                        s = r[0],
                        a = this,
                        o = 1,
                        h = r.length;
                    h > o;
                    o++
                )
                    (i = r[o]), n(s, i), (s = i);
                t._closed && n(i, r[0]);
            },
            _computeParts: function (t, e, n, i) {
                if (i - n > 1 / 32 && !z.isFlatEnough(t, 0.25)) {
                    var r = z.subdivide(t),
                        s = (n + i) / 2;
                    this._computeParts(r[0], e, n, s),
                        this._computeParts(r[1], e, s, i);
                } else {
                    var a = t[6] - t[0],
                        o = t[7] - t[1],
                        h = Math.sqrt(a * a + o * o);
                    h > 1e-5 &&
                        ((this.length += h),
                        this.parts.push({
                            offset: this.length,
                            value: i,
                            index: e
                        }));
                }
            },
            getParameterAt: function (t) {
                for (
                    var e, n = this.index;
                    (e = n), !(0 == n || this.parts[--n].offset < t);

                );
                for (var i = this.parts.length; i > e; e++) {
                    var r = this.parts[e];
                    if (r.offset >= t) {
                        this.index = e;
                        var s = this.parts[e - 1],
                            a = s && s.index == r.index ? s.value : 0,
                            o = s ? s.offset : 0;
                        return {
                            value:
                                a + ((r.value - a) * (t - o)) / (r.offset - o),
                            index: r.index
                        };
                    }
                }
                var r = this.parts[this.parts.length - 1];
                return {
                    value: 1,
                    index: r.index
                };
            },
            evaluate: function (t, e) {
                var n = this.getParameterAt(t);
                return z.evaluate(this.curves[n.index], n.value, e);
            },
            drawPart: function (t, e, n) {
                (e = this.getParameterAt(e)), (n = this.getParameterAt(n));
                for (var i = e.index; i <= n.index; i++) {
                    var r = z.getPart(
                        this.curves[i],
                        i == e.index ? e.value : 0,
                        i == n.index ? n.value : 1
                    );
                    i == e.index && t.moveTo(r[0], r[1]),
                        t.bezierCurveTo.apply(t, r.slice(2));
                }
            }
        }),
        E = e.extend({
            initialize: function (t, e) {
                this.points = [];
                for (var n, i = t._segments, r = 0, s = i.length; s > r; r++) {
                    var a = i[r].point.clone();
                    (n && n.equals(a)) || (this.points.push(a), (n = a));
                }
                this.error = e;
            },
            fit: function () {
                var t = this.points,
                    e = t.length;
                return (
                    (this.segments = e > 0 ? [new k(t[0])] : []),
                    e > 1 &&
                        this.fitCubic(
                            0,
                            e - 1,
                            t[1].subtract(t[0]).normalize(),
                            t[e - 2].subtract(t[e - 1]).normalize()
                        ),
                    this.segments
                );
            },
            fitCubic: function (e, n, i, r) {
                if (1 == n - e) {
                    var s = this.points[e],
                        a = this.points[n],
                        o = s.getDistance(a) / 3;
                    return (
                        this.addCurve([
                            s,
                            s.add(i.normalize(o)),
                            a.add(r.normalize(o)),
                            a
                        ]),
                        t
                    );
                }
                for (
                    var h,
                        u = this.chordLengthParameterize(e, n),
                        l = Math.max(this.error, this.error * this.error),
                        c = 0;
                    4 >= c;
                    c++
                ) {
                    var d = this.generateBezier(e, n, u, i, r),
                        f = this.findMaxError(e, n, d, u);
                    if (f.error < this.error) return this.addCurve(d), t;
                    if (((h = f.index), f.error >= l)) break;
                    this.reparameterize(e, n, u, d), (l = f.error);
                }
                var _ = this.points[h - 1].subtract(this.points[h]),
                    g = this.points[h].subtract(this.points[h + 1]),
                    p = _.add(g).divide(2).normalize();
                this.fitCubic(e, h, i, p), this.fitCubic(h, n, p.negate(), r);
            },
            addCurve: function (t) {
                var e = this.segments[this.segments.length - 1];
                e.setHandleOut(t[1].subtract(t[0])),
                    this.segments.push(new k(t[3], t[2].subtract(t[3])));
            },
            generateBezier: function (t, e, n, i, r) {
                for (
                    var s = 1e-11,
                        a = this.points[t],
                        o = this.points[e],
                        h = [
                            [0, 0],
                            [0, 0]
                        ],
                        u = [0, 0],
                        l = 0,
                        c = e - t + 1;
                    c > l;
                    l++
                ) {
                    var d = n[l],
                        f = 1 - d,
                        _ = 3 * d * f,
                        g = f * f * f,
                        p = _ * f,
                        v = _ * d,
                        m = d * d * d,
                        y = i.normalize(p),
                        w = r.normalize(v),
                        x = this.points[t + l]
                            .subtract(a.multiply(g + p))
                            .subtract(o.multiply(v + m));
                    (h[0][0] += y.dot(y)),
                        (h[0][1] += y.dot(w)),
                        (h[1][0] = h[0][1]),
                        (h[1][1] += w.dot(w)),
                        (u[0] += y.dot(x)),
                        (u[1] += w.dot(x));
                }
                var b,
                    C,
                    S = h[0][0] * h[1][1] - h[1][0] * h[0][1];
                if (Math.abs(S) > s) {
                    var P = h[0][0] * u[1] - h[1][0] * u[0],
                        k = u[0] * h[1][1] - u[1] * h[0][1];
                    (b = k / S), (C = P / S);
                } else {
                    var M = h[0][0] + h[0][1],
                        z = h[1][0] + h[1][1];
                    b = C =
                        Math.abs(M) > s
                            ? u[0] / M
                            : Math.abs(z) > s
                            ? u[1] / z
                            : 0;
                }
                var A = o.getDistance(a);
                return (
                    (s *= A),
                    (s > b || s > C) && (b = C = A / 3),
                    [a, a.add(i.normalize(b)), o.add(r.normalize(C)), o]
                );
            },
            reparameterize: function (t, e, n, i) {
                for (var r = t; e >= r; r++)
                    n[r - t] = this.findRoot(i, this.points[r], n[r - t]);
            },
            findRoot: function (t, e, n) {
                for (var i = [], r = [], s = 0; 2 >= s; s++)
                    i[s] = t[s + 1].subtract(t[s]).multiply(3);
                for (var s = 0; 1 >= s; s++)
                    r[s] = i[s + 1].subtract(i[s]).multiply(2);
                var a = this.evaluate(3, t, n),
                    o = this.evaluate(2, i, n),
                    h = this.evaluate(1, r, n),
                    u = a.subtract(e),
                    l = o.dot(o) + u.dot(h);
                return Math.abs(l) < 1e-5 ? n : n - u.dot(o) / l;
            },
            evaluate: function (t, e, n) {
                for (var i = e.slice(), r = 1; t >= r; r++)
                    for (var s = 0; t - r >= s; s++)
                        i[s] = i[s].multiply(1 - n).add(i[s + 1].multiply(n));
                return i[0];
            },
            chordLengthParameterize: function (t, e) {
                for (var n = [0], i = t + 1; e >= i; i++)
                    n[i - t] =
                        n[i - t - 1] +
                        this.points[i].getDistance(this.points[i - 1]);
                for (var i = 1, r = e - t; r >= i; i++) n[i] /= n[r];
                return n;
            },
            findMaxError: function (t, e, n, i) {
                for (
                    var r = Math.floor((e - t + 1) / 2), s = 0, a = t + 1;
                    e > a;
                    a++
                ) {
                    var o = this.evaluate(3, n, i[a - t]),
                        h = o.subtract(this.points[a]),
                        u = h.x * h.x + h.y * h.y;
                    u >= s && ((s = u), (r = a));
                }
                return {
                    error: s,
                    index: r
                };
            }
        }),
        N = y.extend({
            _class: "TextItem",
            _boundsSelected: !0,
            _applyMatrix: !1,
            _canApplyMatrix: !1,
            _serializeFields: {
                content: null
            },
            _boundsGetter: "getBounds",
            initialize: function (n) {
                (this._content = ""), (this._lines = []);
                var i = n && e.isPlainObject(n) && n.x === t && n.y === t;
                this._initialize(i && n, !i && h.read(arguments));
            },
            _equals: function (t) {
                return this._content === t._content;
            },
            _clone: function be(t) {
                return t.setContent(this._content), be.base.call(this, t);
            },
            getContent: function () {
                return this._content;
            },
            setContent: function (t) {
                (this._content = "" + t),
                    (this._lines = this._content.split(/\r\n|\n|\r/gm)),
                    this._changed(265);
            },
            isEmpty: function () {
                return !this._content;
            },
            getCharacterStyle: "#getStyle",
            setCharacterStyle: "#setStyle",
            getParagraphStyle: "#getStyle",
            setParagraphStyle: "#setStyle"
        }),
        j = N.extend({
            _class: "PointText",
            initialize: function () {
                N.apply(this, arguments);
            },
            clone: function (t) {
                return this._clone(new j(y.NO_INSERT), t);
            },
            getPoint: function () {
                var t = this._matrix.getTranslation();
                return new u(t.x, t.y, this, "setPoint");
            },
            setPoint: function () {
                var t = h.read(arguments);
                this.translate(t.subtract(this._matrix.getTranslation()));
            },
            _draw: function (t) {
                if (this._content) {
                    this._setStyles(t);
                    var e = this._style,
                        n = this._lines,
                        i = e.getLeading(),
                        r = t.shadowColor;
                    (t.font = e.getFontStyle()),
                        (t.textAlign = e.getJustification());
                    for (var s = 0, a = n.length; a > s; s++) {
                        t.shadowColor = r;
                        var o = n[s];
                        e.hasFill() &&
                            (t.fillText(o, 0, 0),
                            (t.shadowColor = "rgba(0,0,0,0)")),
                            e.hasStroke() && t.strokeText(o, 0, 0),
                            t.translate(0, i);
                    }
                }
            },
            _getBounds: function (t, e) {
                var n = this._style,
                    i = this._lines,
                    r = i.length,
                    s = n.getJustification(),
                    a = n.getLeading(),
                    o = this.getView().getTextWidth(n.getFontStyle(), i),
                    h = 0;
                "left" !== s && (h -= o / ("center" === s ? 2 : 1));
                var u = new f(h, r ? -0.75 * a : 0, o, r * a);
                return e ? e._transformBounds(u, u) : u;
            }
        }),
        D = e.extend(
            new (function () {
                function t(t) {
                    var e,
                        i = t.match(/^#(\w{1,2})(\w{1,2})(\w{1,2})$/);
                    if (i) {
                        e = [0, 0, 0];
                        for (var r = 0; 3 > r; r++) {
                            var a = i[r + 1];
                            e[r] =
                                parseInt(1 == a.length ? a + a : a, 16) / 255;
                        }
                    } else if ((i = t.match(/^rgba?\((.*)\)$/))) {
                        e = i[1].split(",");
                        for (var r = 0, o = e.length; o > r; r++) {
                            var a = +e[r];
                            e[r] = 3 > r ? a / 255 : a;
                        }
                    } else {
                        var h = s[t];
                        if (!h) {
                            n ||
                                ((n = Q.getContext(1, 1)),
                                (n.globalCompositeOperation = "copy")),
                                (n.fillStyle = "rgba(0,0,0,0)"),
                                (n.fillStyle = t),
                                n.fillRect(0, 0, 1, 1);
                            var u = n.getImageData(0, 0, 1, 1).data;
                            h = s[t] = [u[0] / 255, u[1] / 255, u[2] / 255];
                        }
                        e = h.slice();
                    }
                    return e;
                }
                var n,
                    i = {
                        gray: ["gray"],
                        rgb: ["red", "green", "blue"],
                        hsb: ["hue", "saturation", "brightness"],
                        hsl: ["hue", "saturation", "lightness"],
                        gradient: [
                            "gradient",
                            "origin",
                            "destination",
                            "highlight"
                        ]
                    },
                    r = {},
                    s = {},
                    o = [
                        [0, 3, 1],
                        [2, 0, 1],
                        [1, 0, 3],
                        [1, 2, 0],
                        [3, 1, 0],
                        [0, 1, 2]
                    ],
                    u = {
                        "rgb-hsb": function (t, e, n) {
                            var i = Math.max(t, e, n),
                                r = Math.min(t, e, n),
                                s = i - r,
                                a =
                                    0 === s
                                        ? 0
                                        : 60 *
                                          (i == t
                                              ? (e - n) / s + (n > e ? 6 : 0)
                                              : i == e
                                              ? (n - t) / s + 2
                                              : (t - e) / s + 4);
                            return [a, 0 === i ? 0 : s / i, i];
                        },
                        "hsb-rgb": function (t, e, n) {
                            t = (((t / 60) % 6) + 6) % 6;
                            var i = Math.floor(t),
                                r = t - i,
                                i = o[i],
                                s = [
                                    n,
                                    n * (1 - e),
                                    n * (1 - e * r),
                                    n * (1 - e * (1 - r))
                                ];
                            return [s[i[0]], s[i[1]], s[i[2]]];
                        },
                        "rgb-hsl": function (t, e, n) {
                            var i = Math.max(t, e, n),
                                r = Math.min(t, e, n),
                                s = i - r,
                                a = 0 === s,
                                o = a
                                    ? 0
                                    : 60 *
                                      (i == t
                                          ? (e - n) / s + (n > e ? 6 : 0)
                                          : i == e
                                          ? (n - t) / s + 2
                                          : (t - e) / s + 4),
                                h = (i + r) / 2,
                                u = a
                                    ? 0
                                    : 0.5 > h
                                    ? s / (i + r)
                                    : s / (2 - i - r);
                            return [o, u, h];
                        },
                        "hsl-rgb": function (t, e, n) {
                            if (((t = (((t / 360) % 1) + 1) % 1), 0 === e))
                                return [n, n, n];
                            for (
                                var i = [t + 1 / 3, t, t - 1 / 3],
                                    r = 0.5 > n ? n * (1 + e) : n + e - n * e,
                                    s = 2 * n - r,
                                    a = [],
                                    o = 0;
                                3 > o;
                                o++
                            ) {
                                var h = i[o];
                                0 > h && (h += 1),
                                    h > 1 && (h -= 1),
                                    (a[o] =
                                        1 > 6 * h
                                            ? s + 6 * (r - s) * h
                                            : 1 > 2 * h
                                            ? r
                                            : 2 > 3 * h
                                            ? s + 6 * (r - s) * (2 / 3 - h)
                                            : s);
                            }
                            return a;
                        },
                        "rgb-gray": function (t, e, n) {
                            return [0.2989 * t + 0.587 * e + 0.114 * n];
                        },
                        "gray-rgb": function (t) {
                            return [t, t, t];
                        },
                        "gray-hsb": function (t) {
                            return [0, 0, t];
                        },
                        "gray-hsl": function (t) {
                            return [0, 0, t];
                        },
                        "gradient-rgb": function () {
                            return [];
                        },
                        "rgb-gradient": function () {
                            return [];
                        }
                    };
                return e.each(
                    i,
                    function (t, n) {
                        (r[n] = []),
                            e.each(
                                t,
                                function (t, s) {
                                    var a = e.capitalize(t),
                                        o = /^(hue|saturation)$/.test(t),
                                        u = (r[n][s] =
                                            "gradient" === t
                                                ? function (t) {
                                                      var e =
                                                          this._components[0];
                                                      return (
                                                          (t = B.read(
                                                              Array.isArray(t)
                                                                  ? t
                                                                  : arguments,
                                                              0,
                                                              {
                                                                  readNull: !0
                                                              }
                                                          )),
                                                          e !== t &&
                                                              (e &&
                                                                  e._removeOwner(
                                                                      this
                                                                  ),
                                                              t &&
                                                                  t._addOwner(
                                                                      this
                                                                  )),
                                                          t
                                                      );
                                                  }
                                                : "gradient" === n
                                                ? function () {
                                                      return h.read(
                                                          arguments,
                                                          0,
                                                          {
                                                              readNull:
                                                                  "highlight" ===
                                                                  t,
                                                              clone: !0
                                                          }
                                                      );
                                                  }
                                                : function (t) {
                                                      return null == t ||
                                                          isNaN(t)
                                                          ? 0
                                                          : t;
                                                  });
                                    (this["get" + a] = function () {
                                        return this._type === n ||
                                            (o && /^hs[bl]$/.test(this._type))
                                            ? this._components[s]
                                            : this._convert(n)[s];
                                    }),
                                        (this["set" + a] = function (t) {
                                            this._type === n ||
                                                (o &&
                                                    /^hs[bl]$/.test(
                                                        this._type
                                                    )) ||
                                                ((this._components =
                                                    this._convert(n)),
                                                (this._properties = i[n]),
                                                (this._type = n)),
                                                (t = u.call(this, t)),
                                                null != t &&
                                                    ((this._components[s] = t),
                                                    this._changed());
                                        });
                                },
                                this
                            );
                    },
                    {
                        _class: "Color",
                        _readIndex: !0,
                        initialize: function l(e) {
                            var n,
                                s,
                                a,
                                o,
                                h = Array.prototype.slice,
                                u = arguments,
                                c = 0;
                            Array.isArray(e) && ((u = e), (e = u[0]));
                            var d = null != e && typeof e;
                            if (
                                ("string" === d &&
                                    e in i &&
                                    ((n = e),
                                    (e = u[1]),
                                    Array.isArray(e)
                                        ? ((s = e), (a = u[2]))
                                        : (this.__read && (c = 1),
                                          (u = h.call(u, 1)),
                                          (d = typeof e))),
                                !s)
                            ) {
                                if (
                                    (o =
                                        "number" === d
                                            ? u
                                            : "object" === d && null != e.length
                                            ? e
                                            : null)
                                ) {
                                    n || (n = o.length >= 3 ? "rgb" : "gray");
                                    var f = i[n].length;
                                    (a = o[f]),
                                        this.__read &&
                                            (c +=
                                                o === arguments
                                                    ? f + (null != a ? 1 : 0)
                                                    : 1),
                                        o.length > f && (o = h.call(o, 0, f));
                                } else if ("string" === d)
                                    (n = "rgb"),
                                        (s = t(e)),
                                        4 === s.length &&
                                            ((a = s[3]), s.length--);
                                else if ("object" === d)
                                    if (e.constructor === l) {
                                        if (
                                            ((n = e._type),
                                            (s = e._components.slice()),
                                            (a = e._alpha),
                                            "gradient" === n)
                                        )
                                            for (
                                                var _ = 1, g = s.length;
                                                g > _;
                                                _++
                                            ) {
                                                var p = s[_];
                                                p && (s[_] = p.clone());
                                            }
                                    } else if (e.constructor === B)
                                        (n = "gradient"), (o = u);
                                    else {
                                        n =
                                            "hue" in e
                                                ? "lightness" in e
                                                    ? "hsl"
                                                    : "hsb"
                                                : "gradient" in e ||
                                                  "stops" in e ||
                                                  "radial" in e
                                                ? "gradient"
                                                : "gray" in e
                                                ? "gray"
                                                : "rgb";
                                        var v = i[n];
                                        (y = r[n]), (this._components = s = []);
                                        for (
                                            var _ = 0, g = v.length;
                                            g > _;
                                            _++
                                        ) {
                                            var m = e[v[_]];
                                            null == m &&
                                                0 === _ &&
                                                "gradient" === n &&
                                                "stops" in e &&
                                                (m = {
                                                    stops: e.stops,
                                                    radial: e.radial
                                                }),
                                                (m = y[_].call(this, m)),
                                                null != m && (s[_] = m);
                                        }
                                        a = e.alpha;
                                    }
                                this.__read && n && (c = 1);
                            }
                            if (
                                ((this._type = n || "rgb"),
                                "gradient" === n &&
                                    (this._id = l._id = (l._id || 0) + 1),
                                !s)
                            ) {
                                this._components = s = [];
                                for (
                                    var y = r[this._type], _ = 0, g = y.length;
                                    g > _;
                                    _++
                                ) {
                                    var m = y[_].call(this, o && o[_]);
                                    null != m && (s[_] = m);
                                }
                            }
                            (this._components = s),
                                (this._properties = i[this._type]),
                                (this._alpha = a),
                                this.__read && (this.__read = c);
                        },
                        _serialize: function (t, n) {
                            var i = this.getComponents();
                            return e.serialize(
                                /^(gray|rgb)$/.test(this._type)
                                    ? i
                                    : [this._type].concat(i),
                                t,
                                !0,
                                n
                            );
                        },
                        _changed: function () {
                            (this._canvasStyle = null),
                                this._owner && this._owner._changed(65);
                        },
                        _convert: function (t) {
                            var e;
                            return this._type === t
                                ? this._components.slice()
                                : (e = u[this._type + "-" + t])
                                ? e.apply(this, this._components)
                                : u["rgb-" + t].apply(
                                      this,
                                      u[this._type + "-rgb"].apply(
                                          this,
                                          this._components
                                      )
                                  );
                        },
                        convert: function (t) {
                            return new D(t, this._convert(t), this._alpha);
                        },
                        getType: function () {
                            return this._type;
                        },
                        setType: function (t) {
                            (this._components = this._convert(t)),
                                (this._properties = i[t]),
                                (this._type = t);
                        },
                        getComponents: function () {
                            var t = this._components.slice();
                            return (
                                null != this._alpha && t.push(this._alpha), t
                            );
                        },
                        getAlpha: function () {
                            return null != this._alpha ? this._alpha : 1;
                        },
                        setAlpha: function (t) {
                            (this._alpha =
                                null == t ? null : Math.min(Math.max(t, 0), 1)),
                                this._changed();
                        },
                        hasAlpha: function () {
                            return null != this._alpha;
                        },
                        equals: function (t) {
                            var n = e.isPlainValue(t, !0)
                                ? D.read(arguments)
                                : t;
                            return (
                                n === this ||
                                (n &&
                                    this._class === n._class &&
                                    this._type === n._type &&
                                    this._alpha === n._alpha &&
                                    e.equals(
                                        this._components,
                                        n._components
                                    )) ||
                                !1
                            );
                        },
                        toString: function () {
                            for (
                                var t = this._properties,
                                    e = [],
                                    n = "gradient" === this._type,
                                    i = a.instance,
                                    r = 0,
                                    s = t.length;
                                s > r;
                                r++
                            ) {
                                var o = this._components[r];
                                null != o &&
                                    e.push(t[r] + ": " + (n ? o : i.number(o)));
                            }
                            return (
                                null != this._alpha &&
                                    e.push("alpha: " + i.number(this._alpha)),
                                "{ " + e.join(", ") + " }"
                            );
                        },
                        toCSS: function (t) {
                            function e(t) {
                                return Math.round(
                                    255 * (0 > t ? 0 : t > 1 ? 1 : t)
                                );
                            }
                            var n = this._convert("rgb"),
                                i = t || null == this._alpha ? 1 : this._alpha;
                            return (
                                (n = [e(n[0]), e(n[1]), e(n[2])]),
                                1 > i && n.push(0 > i ? 0 : i),
                                t
                                    ? "#" +
                                      (
                                          (1 << 24) +
                                          (n[0] << 16) +
                                          (n[1] << 8) +
                                          n[2]
                                      )
                                          .toString(16)
                                          .slice(1)
                                    : (4 == n.length ? "rgba(" : "rgb(") +
                                      n.join(",") +
                                      ")"
                            );
                        },
                        toCanvasStyle: function (t) {
                            if (this._canvasStyle) return this._canvasStyle;
                            if ("gradient" !== this._type)
                                return (this._canvasStyle = this.toCSS());
                            var e,
                                n = this._components,
                                i = n[0],
                                r = i._stops,
                                s = n[1],
                                a = n[2];
                            if (i._radial) {
                                var o = a.getDistance(s),
                                    h = n[3];
                                if (h) {
                                    var u = h.subtract(s);
                                    u.getLength() > o &&
                                        (h = s.add(u.normalize(o - 0.1)));
                                }
                                var l = h || s;
                                e = t.createRadialGradient(
                                    l.x,
                                    l.y,
                                    0,
                                    s.x,
                                    s.y,
                                    o
                                );
                            } else
                                e = t.createLinearGradient(s.x, s.y, a.x, a.y);
                            for (var c = 0, d = r.length; d > c; c++) {
                                var f = r[c];
                                e.addColorStop(
                                    f._rampPoint,
                                    f._color.toCanvasStyle()
                                );
                            }
                            return (this._canvasStyle = e);
                        },
                        transform: function (t) {
                            if ("gradient" === this._type) {
                                for (
                                    var e = this._components,
                                        n = 1,
                                        i = e.length;
                                    i > n;
                                    n++
                                ) {
                                    var r = e[n];
                                    t._transformPoint(r, r, !0);
                                }
                                this._changed();
                            }
                        },
                        statics: {
                            _types: i,
                            random: function () {
                                var t = Math.random;
                                return new D(t(), t(), t());
                            }
                        }
                    }
                );
            })(),
            new (function () {
                var t = {
                    add: function (t, e) {
                        return t + e;
                    },
                    subtract: function (t, e) {
                        return t - e;
                    },
                    multiply: function (t, e) {
                        return t * e;
                    },
                    divide: function (t, e) {
                        return t / e;
                    }
                };
                return e.each(
                    t,
                    function (t, e) {
                        this[e] = function (e) {
                            e = D.read(arguments);
                            for (
                                var n = this._type,
                                    i = this._components,
                                    r = e._convert(n),
                                    s = 0,
                                    a = i.length;
                                a > s;
                                s++
                            )
                                r[s] = t(i[s], r[s]);
                            return new D(
                                n,
                                r,
                                null != this._alpha
                                    ? t(this._alpha, e.getAlpha())
                                    : null
                            );
                        };
                    },
                    {}
                );
            })()
        );
    e.each(
        D._types,
        function (t, n) {
            var i = (this[e.capitalize(n) + "Color"] = function (t) {
                var e = null != t && typeof t,
                    i =
                        "object" === e && null != t.length
                            ? t
                            : "string" === e
                            ? null
                            : arguments;
                return i ? new D(n, i) : new D(t);
            });
            if (3 == n.length) {
                var r = n.toUpperCase();
                D[r] = this[r + "Color"] = i;
            }
        },
        e.exports
    );
    var B = e.extend({
            _class: "Gradient",
            initialize: function Ce(t, e) {
                (this._id = Ce._id = (Ce._id || 0) + 1),
                    t && this._set(t) && (t = e = null),
                    this._stops || this.setStops(t || ["white", "black"]),
                    null == this._radial &&
                        this.setRadial(
                            ("string" == typeof e && "radial" === e) || e || !1
                        );
            },
            _serialize: function (t, n) {
                return n.add(this, function () {
                    return e.serialize([this._stops, this._radial], t, !0, n);
                });
            },
            _changed: function () {
                for (
                    var t = 0, e = this._owners && this._owners.length;
                    e > t;
                    t++
                )
                    this._owners[t]._changed();
            },
            _addOwner: function (t) {
                this._owners || (this._owners = []), this._owners.push(t);
            },
            _removeOwner: function (e) {
                var n = this._owners ? this._owners.indexOf(e) : -1;
                -1 != n &&
                    (this._owners.splice(n, 1),
                    0 === this._owners.length && (this._owners = t));
            },
            clone: function () {
                for (var t = [], e = 0, n = this._stops.length; n > e; e++)
                    t[e] = this._stops[e].clone();
                return new B(t);
            },
            getStops: function () {
                return this._stops;
            },
            setStops: function (e) {
                if (this.stops)
                    for (var n = 0, i = this._stops.length; i > n; n++)
                        this._stops[n]._owner = t;
                if (e.length < 2)
                    throw Error(
                        "Gradient stop list needs to contain at least two stops."
                    );
                this._stops = R.readAll(e, 0, {
                    clone: !0
                });
                for (var n = 0, i = this._stops.length; i > n; n++) {
                    var r = this._stops[n];
                    (r._owner = this),
                        r._defaultRamp && r.setRampPoint(n / (i - 1));
                }
                this._changed();
            },
            getRadial: function () {
                return this._radial;
            },
            setRadial: function (t) {
                (this._radial = t), this._changed();
            },
            equals: function (t) {
                if (t === this) return !0;
                if (
                    t &&
                    this._class === t._class &&
                    this._stops.length === t._stops.length
                ) {
                    for (var e = 0, n = this._stops.length; n > e; e++)
                        if (!this._stops[e].equals(t._stops[e])) return !1;
                    return !0;
                }
                return !1;
            }
        }),
        R = e.extend({
            _class: "GradientStop",
            initialize: function (e, n) {
                if (e) {
                    var i, r;
                    n === t && Array.isArray(e)
                        ? ((i = e[0]), (r = e[1]))
                        : e.color
                        ? ((i = e.color), (r = e.rampPoint))
                        : ((i = e), (r = n)),
                        this.setColor(i),
                        this.setRampPoint(r);
                }
            },
            clone: function () {
                return new R(this._color.clone(), this._rampPoint);
            },
            _serialize: function (t, n) {
                return e.serialize([this._color, this._rampPoint], t, !0, n);
            },
            _changed: function () {
                this._owner && this._owner._changed(65);
            },
            getRampPoint: function () {
                return this._rampPoint;
            },
            setRampPoint: function (t) {
                (this._defaultRamp = null == t),
                    (this._rampPoint = t || 0),
                    this._changed();
            },
            getColor: function () {
                return this._color;
            },
            setColor: function (t) {
                (this._color = D.read(arguments)),
                    this._color === t && (this._color = t.clone()),
                    (this._color._owner = this),
                    this._changed();
            },
            equals: function (t) {
                return (
                    t === this ||
                    (t &&
                        this._class === t._class &&
                        this._color.equals(t._color) &&
                        this._rampPoint == t._rampPoint) ||
                    !1
                );
            }
        }),
        F = e.extend(
            new (function () {
                var n = {
                        fillColor: t,
                        strokeColor: t,
                        strokeWidth: 1,
                        strokeCap: "butt",
                        strokeJoin: "miter",
                        strokeScaling: !0,
                        miterLimit: 10,
                        dashOffset: 0,
                        dashArray: [],
                        windingRule: "nonzero",
                        shadowColor: t,
                        shadowBlur: 0,
                        shadowOffset: new h(),
                        selectedColor: t,
                        fontFamily: "sans-serif",
                        fontWeight: "normal",
                        fontSize: 12,
                        font: "sans-serif",
                        leading: null,
                        justification: "left"
                    },
                    i = {
                        strokeWidth: 97,
                        strokeCap: 97,
                        strokeJoin: 97,
                        strokeScaling: 105,
                        miterLimit: 97,
                        fontFamily: 9,
                        fontWeight: 9,
                        fontSize: 9,
                        font: 9,
                        leading: 9,
                        justification: 9
                    },
                    r = {
                        beans: !0
                    },
                    s = {
                        _defaults: n,
                        _textDefaults: new e(n, {
                            fillColor: new D()
                        }),
                        beans: !0
                    };
                return (
                    e.each(n, function (n, a) {
                        var o = /Color$/.test(a),
                            h = e.capitalize(a),
                            u = i[a],
                            l = "set" + h,
                            c = "get" + h;
                        (s[l] = function (e) {
                            var n = this._owner,
                                i = n && n._children;
                            if (i && i.length > 0 && !(n instanceof T))
                                for (var r = 0, s = i.length; s > r; r++)
                                    i[r]._style[l](e);
                            else {
                                var h = this._values[a];
                                h != e &&
                                    (o &&
                                        (h && (h._owner = t),
                                        e &&
                                            e.constructor === D &&
                                            (e._owner && (e = e.clone()),
                                            (e._owner = n))),
                                    (this._values[a] = e),
                                    n && n._changed(u || 65));
                            }
                        }),
                            (s[c] = function (n) {
                                var i,
                                    r = this._owner,
                                    s = r && r._children;
                                if (
                                    !s ||
                                    0 === s.length ||
                                    n ||
                                    r instanceof T
                                ) {
                                    var i = this._values[a];
                                    return (
                                        i === t
                                            ? ((i = this._defaults[a]),
                                              i && i.clone && (i = i.clone()),
                                              (this._values[a] = i))
                                            : !o ||
                                              (i && i.constructor === D) ||
                                              ((this._values[a] = i =
                                                  D.read([i], 0, {
                                                      readNull: !0,
                                                      clone: !0
                                                  })),
                                              i && (i._owner = r)),
                                        i
                                    );
                                }
                                for (var h = 0, u = s.length; u > h; h++) {
                                    var l = s[h]._style[c]();
                                    if (0 === h) i = l;
                                    else if (!e.equals(i, l)) return t;
                                }
                                return i;
                            }),
                            (r[c] = function (t) {
                                return this._style[c](t);
                            }),
                            (r[l] = function (t) {
                                this._style[l](t);
                            });
                    }),
                    y.inject(r),
                    s
                );
            })(),
            {
                _class: "Style",
                initialize: function (t, e, n) {
                    (this._values = {}),
                        (this._owner = e),
                        (this._project =
                            (e && e._project) || n || paper.project),
                        e instanceof N && (this._defaults = this._textDefaults),
                        t && this.set(t);
                },
                set: function (t) {
                    var e = t instanceof F,
                        n = e ? t._values : t;
                    if (n)
                        for (var i in n)
                            if (i in this._defaults) {
                                var r = n[i];
                                this[i] = r && e && r.clone ? r.clone() : r;
                            }
                },
                equals: function (t) {
                    return (
                        t === this ||
                        (t &&
                            this._class === t._class &&
                            e.equals(this._values, t._values)) ||
                        !1
                    );
                },
                hasFill: function () {
                    return !!this.getFillColor();
                },
                hasStroke: function () {
                    return !!this.getStrokeColor() && this.getStrokeWidth() > 0;
                },
                hasShadow: function () {
                    return !!this.getShadowColor() && this.getShadowBlur() > 0;
                },
                getView: function () {
                    return this._project.getView();
                },
                getFontStyle: function () {
                    var t = this.getFontSize();
                    return (
                        this.getFontWeight() +
                        " " +
                        t +
                        (/[a-z]/i.test(t + "") ? " " : "px ") +
                        this.getFontFamily()
                    );
                },
                getFont: "#getFontFamily",
                setFont: "#setFontFamily",
                getLeading: function Se() {
                    var t = Se.base.call(this),
                        e = this.getFontSize();
                    return (
                        /pt|em|%|px/.test(e) &&
                            (e = this.getView().getPixelSize(e)),
                        null != t ? t : 1.2 * e
                    );
                }
            }
        ),
        V = new (function () {
            function n(t, i) {
                for (var r = [], s = 0, a = t && t.length; a > s; ) {
                    var o = t[s++];
                    if ("string" == typeof o) o = document.createElement(o);
                    else if (!o || !o.nodeType) continue;
                    e.isPlainObject(t[s]) && V.set(o, t[s++]),
                        Array.isArray(t[s]) && n(t[s++], o),
                        i && i.appendChild(o),
                        r.push(o);
                }
                return r;
            }
            function i(t, e, n, i) {
                for (
                    var r = ["", "webkit", "moz", "Moz", "ms", "o"],
                        s = e[0].toUpperCase() + e.substring(1),
                        a = 0;
                    6 > a;
                    a++
                ) {
                    var o = r[a],
                        h = o ? o + s : e;
                    if (h in t) {
                        if (!n) return t[h];
                        t[h] = i;
                        break;
                    }
                }
            }
            var r = /^(checked|value|selected|disabled)$/i,
                s = {
                    text: "textContent",
                    html: "innerHTML"
                },
                a = {
                    lineHeight: 1,
                    zoom: 1,
                    zIndex: 1,
                    opacity: 1
                };
            return {
                create: function (t, e) {
                    var i = Array.isArray(t),
                        r = n(i ? t : arguments, i ? e : null);
                    return 1 == r.length ? r[0] : r;
                },
                find: function (t, e) {
                    return (e || document).querySelector(t);
                },
                findAll: function (t, e) {
                    return (e || document).querySelectorAll(t);
                },
                get: function (t, e) {
                    return t
                        ? r.test(e)
                            ? "value" === e || "string" != typeof t[e]
                                ? t[e]
                                : !0
                            : e in s
                            ? t[s[e]]
                            : t.getAttribute(e)
                        : null;
                },
                set: function (e, n, i) {
                    if ("string" != typeof n)
                        for (var a in n)
                            n.hasOwnProperty(a) && this.set(e, a, n[a]);
                    else {
                        if (!e || i === t) return e;
                        r.test(n)
                            ? (e[n] = i)
                            : n in s
                            ? (e[s[n]] = i)
                            : "style" === n
                            ? this.setStyle(e, i)
                            : "events" === n
                            ? q.add(e, i)
                            : e.setAttribute(n, i);
                    }
                    return e;
                },
                getStyles: function (t) {
                    var e = t && 9 !== t.nodeType ? t.ownerDocument : t,
                        n = e && e.defaultView;
                    return n && n.getComputedStyle(t, "");
                },
                getStyle: function (t, e) {
                    return (t && t.style[e]) || this.getStyles(t)[e] || null;
                },
                setStyle: function (t, e, n) {
                    if ("string" != typeof e)
                        for (var i in e)
                            e.hasOwnProperty(i) && this.setStyle(t, i, e[i]);
                    else
                        !/^-?[\d\.]+$/.test(n) || e in a || (n += "px"),
                            (t.style[e] = n);
                    return t;
                },
                hasClass: function (t, e) {
                    return RegExp("\\s*" + e + "\\s*").test(t.className);
                },
                addClass: function (t, e) {
                    t.className = (t.className + " " + e).trim();
                },
                removeClass: function (t, e) {
                    t.className = t.className
                        .replace(RegExp("\\s*" + e + "\\s*"), " ")
                        .trim();
                },
                remove: function (t) {
                    t.parentNode && t.parentNode.removeChild(t);
                },
                removeChildren: function (t) {
                    for (; t.firstChild; ) t.removeChild(t.firstChild);
                },
                getBounds: function (t, e) {
                    var n,
                        i = t.ownerDocument,
                        r = i.body,
                        s = i.documentElement;
                    try {
                        n = t.getBoundingClientRect();
                    } catch (a) {
                        n = {
                            left: 0,
                            top: 0,
                            width: 0,
                            height: 0
                        };
                    }
                    var o = n.left - (s.clientLeft || r.clientLeft || 0),
                        h = n.top - (s.clientTop || r.clientTop || 0);
                    if (!e) {
                        var u = i.defaultView;
                        (o += u.pageXOffset || s.scrollLeft || r.scrollLeft),
                            (h += u.pageYOffset || s.scrollTop || r.scrollTop);
                    }
                    return new f(o, h, n.width, n.height);
                },
                getViewportBounds: function (t) {
                    var e = t.ownerDocument,
                        n = e.defaultView,
                        i = e.documentElement;
                    return new f(
                        0,
                        0,
                        n.innerWidth || i.clientWidth,
                        n.innerHeight || i.clientHeight
                    );
                },
                getOffset: function (t, e) {
                    return this.getBounds(t, e).getPoint();
                },
                getSize: function (t) {
                    return this.getBounds(t, !0).getSize();
                },
                isInvisible: function (t) {
                    return this.getSize(t).equals(new c(0, 0));
                },
                isInView: function (t) {
                    return (
                        !this.isInvisible(t) &&
                        this.getViewportBounds(t).intersects(
                            this.getBounds(t, !0)
                        )
                    );
                },
                getPrefixed: function (t, e) {
                    return i(t, e);
                },
                setPrefixed: function (t, e, n) {
                    if ("object" == typeof e)
                        for (var r in e) i(t, r, !0, e[r]);
                    else i(t, e, !0, n);
                }
            };
        })(),
        q = {
            add: function (t, e) {
                for (var n in e)
                    for (
                        var i = e[n],
                            r = n.split(/[\s,]+/g),
                            s = 0,
                            a = r.length;
                        a > s;
                        s++
                    )
                        t.addEventListener(r[s], i, !1);
            },
            remove: function (t, e) {
                for (var n in e)
                    for (
                        var i = e[n],
                            r = n.split(/[\s,]+/g),
                            s = 0,
                            a = r.length;
                        a > s;
                        s++
                    )
                        t.removeEventListener(r[s], i, !1);
            },
            getPoint: function (t) {
                var e = t.targetTouches
                    ? t.targetTouches.length
                        ? t.targetTouches[0]
                        : t.changedTouches[0]
                    : t;
                return new h(
                    e.pageX || e.clientX + document.documentElement.scrollLeft,
                    e.pageY || e.clientY + document.documentElement.scrollTop
                );
            },
            getTarget: function (t) {
                return t.target || t.srcElement;
            },
            getRelatedTarget: function (t) {
                return t.relatedTarget || t.toElement;
            },
            getOffset: function (t, e) {
                return q.getPoint(t).subtract(V.getOffset(e || q.getTarget(t)));
            },
            stop: function (t) {
                t.stopPropagation(), t.preventDefault();
            }
        };
    q.requestAnimationFrame = new (function () {
        function t() {
            for (var e = s.length - 1; e >= 0; e--) {
                var o = s[e],
                    h = o[0],
                    u = o[1];
                (!u ||
                    (("true" == r.getAttribute(u, "keepalive") || a) &&
                        V.isInView(u))) &&
                    (s.splice(e, 1), h());
            }
            n && (s.length ? n(t) : (i = !1));
        }
        var e,
            n = V.getPrefixed(window, "requestAnimationFrame"),
            i = !1,
            s = [],
            a = !0;
        return (
            q.add(window, {
                focus: function () {
                    a = !0;
                },
                blur: function () {
                    a = !1;
                }
            }),
            function (r, a) {
                s.push([r, a]),
                    n
                        ? i || (n(t), (i = !0))
                        : e || (e = setInterval(t, 1e3 / 60));
            }
        );
    })();
    var Z = e.extend(
            n,
            {
                _class: "View",
                initialize: function Pe(t, e) {
                    (this._project = t),
                        (this._scope = t._scope),
                        (this._element = e);
                    var n;
                    this._pixelRatio ||
                        (this._pixelRatio = window.devicePixelRatio || 1),
                        (this._id = e.getAttribute("id")),
                        null == this._id &&
                            e.setAttribute(
                                "id",
                                (this._id = "view-" + Pe._id++)
                            ),
                        q.add(e, this._viewEvents);
                    var i = "none";
                    if (
                        (V.setPrefixed(e.style, {
                            userSelect: i,
                            touchAction: i,
                            touchCallout: i,
                            contentZooming: i,
                            userDrag: i,
                            tapHighlightColor: "rgba(0,0,0,0)"
                        }),
                        r.hasAttribute(e, "resize"))
                    ) {
                        var s = V.getOffset(e, !0),
                            a = this;
                        (n = V.getViewportBounds(e).getSize().subtract(s)),
                            (this._windowEvents = {
                                resize: function () {
                                    V.isInvisible(e) ||
                                        (s = V.getOffset(e, !0)),
                                        a.setViewSize(
                                            V.getViewportBounds(e)
                                                .getSize()
                                                .subtract(s)
                                        );
                                }
                            }),
                            q.add(window, this._windowEvents);
                    } else if (((n = V.getSize(e)), n.isNaN() || n.isZero())) {
                        var o = function (t) {
                            return e[t] || parseInt(e.getAttribute(t), 10);
                        };
                        n = new c(o("width"), o("height"));
                    }
                    if (
                        (this._setViewSize(n),
                        r.hasAttribute(e, "stats") &&
                            "undefined" != typeof Stats)
                    ) {
                        this._stats = new Stats();
                        var h = this._stats.domElement,
                            u = h.style,
                            s = V.getOffset(e);
                        (u.position = "absolute"),
                            (u.left = s.x + "px"),
                            (u.top = s.y + "px"),
                            document.body.appendChild(h);
                    }
                    Pe._views.push(this),
                        (Pe._viewsById[this._id] = this),
                        (this._viewSize = n),
                        ((this._matrix = new g())._owner = this),
                        (this._zoom = 1),
                        Pe._focused || (Pe._focused = this),
                        (this._frameItems = {}),
                        (this._frameItemCount = 0);
                },
                remove: function () {
                    return this._project
                        ? (Z._focused === this && (Z._focused = null),
                          Z._views.splice(Z._views.indexOf(this), 1),
                          delete Z._viewsById[this._id],
                          this._project._view === this &&
                              (this._project._view = null),
                          q.remove(this._element, this._viewEvents),
                          q.remove(window, this._windowEvents),
                          (this._element = this._project = null),
                          this.detach("frame"),
                          (this._animate = !1),
                          (this._frameItems = {}),
                          !0)
                        : !1;
                },
                _events: {
                    onFrame: {
                        install: function () {
                            this.play();
                        },
                        uninstall: function () {
                            this.pause();
                        }
                    },
                    onResize: {}
                },
                _animate: !1,
                _time: 0,
                _count: 0,
                _requestFrame: function () {
                    var t = this;
                    q.requestAnimationFrame(function () {
                        (t._requested = !1),
                            t._animate && (t._requestFrame(), t._handleFrame());
                    }, this._element),
                        (this._requested = !0);
                },
                _handleFrame: function () {
                    paper = this._scope;
                    var t = Date.now() / 1e3,
                        n = this._before ? t - this._before : 0;
                    (this._before = t),
                        (this._handlingFrame = !0),
                        this.fire(
                            "frame",
                            new e({
                                delta: n,
                                time: (this._time += n),
                                count: this._count++
                            })
                        ),
                        this._stats && this._stats.update(),
                        (this._handlingFrame = !1),
                        this.update();
                },
                _animateItem: function (t, e) {
                    var n = this._frameItems;
                    e
                        ? ((n[t._id] = {
                              item: t,
                              time: 0,
                              count: 0
                          }),
                          1 === ++this._frameItemCount &&
                              this.attach("frame", this._handleFrameItems))
                        : (delete n[t._id],
                          0 === --this._frameItemCount &&
                              this.detach("frame", this._handleFrameItems));
                },
                _handleFrameItems: function (t) {
                    for (var n in this._frameItems) {
                        var i = this._frameItems[n];
                        i.item.fire(
                            "frame",
                            new e(t, {
                                time: (i.time += t.delta),
                                count: i.count++
                            })
                        );
                    }
                },
                _update: function () {
                    (this._project._needsUpdate = !0),
                        this._handlingFrame ||
                            (this._animate
                                ? this._handleFrame()
                                : this.update());
                },
                _changed: function (t) {
                    1 & t && (this._project._needsUpdate = !0);
                },
                _transform: function (t) {
                    this._matrix.concatenate(t),
                        (this._bounds = null),
                        this._update();
                },
                getElement: function () {
                    return this._element;
                },
                getPixelRatio: function () {
                    return this._pixelRatio;
                },
                getResolution: function () {
                    return 72 * this._pixelRatio;
                },
                getViewSize: function () {
                    var t = this._viewSize;
                    return new d(t.width, t.height, this, "setViewSize");
                },
                setViewSize: function () {
                    var t = c.read(arguments),
                        e = t.subtract(this._viewSize);
                    e.isZero() ||
                        (this._viewSize.set(t.width, t.height),
                        this._setViewSize(t),
                        (this._bounds = null),
                        this.fire("resize", {
                            size: t,
                            delta: e
                        }),
                        this._update());
                },
                _setViewSize: function (t) {
                    var e = this._element;
                    (e.width = t.width), (e.height = t.height);
                },
                getBounds: function () {
                    return (
                        this._bounds ||
                            (this._bounds = this._matrix
                                .inverted()
                                ._transformBounds(
                                    new f(new h(), this._viewSize)
                                )),
                        this._bounds
                    );
                },
                getSize: function () {
                    return this.getBounds().getSize();
                },
                getCenter: function () {
                    return this.getBounds().getCenter();
                },
                setCenter: function () {
                    var t = h.read(arguments);
                    this.scrollBy(t.subtract(this.getCenter()));
                },
                getZoom: function () {
                    return this._zoom;
                },
                setZoom: function (t) {
                    this._transform(
                        new g().scale(t / this._zoom, this.getCenter())
                    ),
                        (this._zoom = t);
                },
                isVisible: function () {
                    return V.isInView(this._element);
                },
                scrollBy: function () {
                    this._transform(
                        new g().translate(h.read(arguments).negate())
                    );
                },
                play: function () {
                    (this._animate = !0),
                        this._requested || this._requestFrame();
                },
                pause: function () {
                    this._animate = !1;
                },
                draw: function () {
                    this.update();
                },
                projectToView: function () {
                    return this._matrix._transformPoint(h.read(arguments));
                },
                viewToProject: function () {
                    return this._matrix._inverseTransform(h.read(arguments));
                }
            },
            {
                statics: {
                    _views: [],
                    _viewsById: {},
                    _id: 0,
                    create: function (t, e) {
                        return (
                            "string" == typeof e &&
                                (e = document.getElementById(e)),
                            new U(t, e)
                        );
                    }
                }
            },
            new (function () {
                function t(t) {
                    var e = q.getTarget(t);
                    return e.getAttribute && Z._viewsById[e.getAttribute("id")];
                }
                function e(t, e) {
                    return t.viewToProject(q.getOffset(e, t._element));
                }
                function n() {
                    if (!Z._focused || !Z._focused.isVisible())
                        for (var t = 0, e = Z._views.length; e > t; t++) {
                            var n = Z._views[t];
                            if (n && n.isVisible()) {
                                Z._focused = a = n;
                                break;
                            }
                        }
                }
                function i(t, e, n) {
                    t._handleEvent("mousemove", e, n);
                    var i = t._scope.tool;
                    return (
                        i &&
                            i._handleEvent(
                                l && i.responds("mousedrag")
                                    ? "mousedrag"
                                    : "mousemove",
                                e,
                                n
                            ),
                        t.update(),
                        i
                    );
                }
                var r,
                    s,
                    a,
                    o,
                    h,
                    u,
                    l = !1,
                    c = window.navigator;
                c.pointerEnabled || c.msPointerEnabled
                    ? ((o = "pointerdown MSPointerDown"),
                      (h = "pointermove MSPointerMove"),
                      (u =
                          "pointerup pointercancel MSPointerUp MSPointerCancel"))
                    : ((o = "touchstart"),
                      (h = "touchmove"),
                      (u = "touchend touchcancel"),
                      ("ontouchstart" in window &&
                          c.userAgent.match(
                              /mobile|tablet|ip(ad|hone|od)|android|silk/i
                          )) ||
                          ((o += " mousedown"),
                          (h += " mousemove"),
                          (u += " mouseup")));
                var d = {
                        "selectstart dragstart": function (t) {
                            l && t.preventDefault();
                        }
                    },
                    f = {
                        mouseout: function (t) {
                            var n = Z._focused,
                                r = q.getRelatedTarget(t);
                            !n ||
                                (r && "HTML" !== r.nodeName) ||
                                i(n, e(n, t), t);
                        },
                        scroll: n
                    };
                return (
                    (d[o] = function (n) {
                        var i = (Z._focused = t(n)),
                            s = e(i, n);
                        (l = !0),
                            i._handleEvent("mousedown", s, n),
                            (r = i._scope.tool) &&
                                r._handleEvent("mousedown", s, n),
                            i.update();
                    }),
                    (f[h] = function (o) {
                        var h = Z._focused;
                        if (!l) {
                            var u = t(o);
                            u
                                ? (h !== u && i(h, e(h, o), o),
                                  (s = h),
                                  (h = Z._focused = a = u))
                                : a && a === h && ((h = Z._focused = s), n());
                        }
                        if (h) {
                            var c = e(h, o);
                            (l || h.getBounds().contains(c)) &&
                                (r = i(h, c, o));
                        }
                    }),
                    (f[u] = function (t) {
                        var n = Z._focused;
                        if (n && l) {
                            var i = e(n, t);
                            (l = !1),
                                n._handleEvent("mouseup", i, t),
                                r && r._handleEvent("mouseup", i, t),
                                n.update();
                        }
                    }),
                    q.add(document, f),
                    q.add(window, {
                        load: n
                    }),
                    {
                        _viewEvents: d,
                        _handleEvent: function () {},
                        statics: {
                            updateFocus: n
                        }
                    }
                );
            })()
        ),
        U = Z.extend(
            {
                _class: "CanvasView",
                initialize: function (t, e) {
                    if (!(e instanceof HTMLCanvasElement)) {
                        var n = c.read(arguments);
                        if (n.isZero())
                            throw Error(
                                "Cannot create CanvasView with the provided argument: " +
                                    [].slice.call(arguments, 1)
                            );
                        e = Q.getCanvas(n);
                    }
                    if (
                        ((this._context = e.getContext("2d")),
                        (this._eventCounters = {}),
                        (this._pixelRatio = 1),
                        "off" !== r.getAttribute(e, "hidpi"))
                    ) {
                        var i = window.devicePixelRatio || 1,
                            s =
                                V.getPrefixed(
                                    this._context,
                                    "backingStorePixelRatio"
                                ) || 1;
                        this._pixelRatio = i / s;
                    }
                    Z.call(this, t, e);
                },
                _setViewSize: function (t) {
                    var e = t.width,
                        n = t.height,
                        i = this._pixelRatio,
                        r = this._element,
                        s = r.style;
                    (r.width = e * i),
                        (r.height = n * i),
                        1 !== i &&
                            ((s.width = e + "px"),
                            (s.height = n + "px"),
                            this._context.scale(i, i));
                },
                getPixelSize: function (t) {
                    var e = this._context,
                        n = e.font;
                    return (
                        (e.font = t + " serif"),
                        (t = parseFloat(e.font)),
                        (e.font = n),
                        t
                    );
                },
                getTextWidth: function (t, e) {
                    var n = this._context,
                        i = n.font,
                        r = 0;
                    n.font = t;
                    for (var s = 0, a = e.length; a > s; s++)
                        r = Math.max(r, n.measureText(e[s]).width);
                    return (n.font = i), r;
                },
                update: function () {
                    var t = this._project;
                    if (!t || !t._needsUpdate) return !1;
                    var e = this._context,
                        n = this._viewSize;
                    return (
                        e.clearRect(0, 0, n.width + 1, n.height + 1),
                        t.draw(e, this._matrix, this._pixelRatio),
                        (t._needsUpdate = !1),
                        !0
                    );
                }
            },
            new (function () {
                function e(e, n, i, r, s, a) {
                    function o(e) {
                        return e.responds(n) &&
                            (h ||
                                (h = new $(
                                    n,
                                    i,
                                    r,
                                    s,
                                    a ? r.subtract(a) : null
                                )),
                            e.fire(n, h) && h.isStopped)
                            ? (i.preventDefault(), !0)
                            : t;
                    }
                    for (var h, u = s; u; ) {
                        if (o(u)) return !0;
                        u = u.getParent();
                    }
                    return o(e) ? !0 : !1;
                }
                var n, i, r, s, a, o, h, u, l;
                return {
                    _handleEvent: function (t, c, d) {
                        if (this._eventCounters[t]) {
                            var f = this._project,
                                _ = f.hitTest(c, {
                                    tolerance: 0,
                                    fill: !0,
                                    stroke: !0
                                }),
                                g = _ && _.item,
                                p = !1;
                            switch (t) {
                                case "mousedown":
                                    for (
                                        p = e(this, t, d, c, g),
                                            u = a == g && Date.now() - l < 300,
                                            s = a = g,
                                            n = i = r = c,
                                            h = !p && g;
                                        h && !h.responds("mousedrag");

                                    )
                                        h = h._parent;
                                    break;
                                case "mouseup":
                                    (p = e(this, t, d, c, g, n)),
                                        h &&
                                            (i &&
                                                !i.equals(c) &&
                                                e(
                                                    this,
                                                    "mousedrag",
                                                    d,
                                                    c,
                                                    h,
                                                    i
                                                ),
                                            g !== h &&
                                                ((r = c),
                                                e(
                                                    this,
                                                    "mousemove",
                                                    d,
                                                    c,
                                                    g,
                                                    r
                                                ))),
                                        !p &&
                                            g &&
                                            g === s &&
                                            ((l = Date.now()),
                                            e(
                                                this,
                                                u && s.responds("doubleclick")
                                                    ? "doubleclick"
                                                    : "click",
                                                d,
                                                n,
                                                g
                                            ),
                                            (u = !1)),
                                        (s = h = null);
                                    break;
                                case "mousemove":
                                    h && (p = e(this, "mousedrag", d, c, h, i)),
                                        p ||
                                            (g !== o && (r = c),
                                            (p = e(this, t, d, c, g, r))),
                                        (i = r = c),
                                        g !== o &&
                                            (e(this, "mouseleave", d, c, o),
                                            (o = g),
                                            e(this, "mouseenter", d, c, g));
                            }
                            return p;
                        }
                    }
                };
            })()
        ),
        H = e.extend({
            _class: "Event",
            initialize: function (t) {
                this.event = t;
            },
            isPrevented: !1,
            isStopped: !1,
            preventDefault: function () {
                (this.isPrevented = !0), this.event.preventDefault();
            },
            stopPropagation: function () {
                (this.isStopped = !0), this.event.stopPropagation();
            },
            stop: function () {
                this.stopPropagation(), this.preventDefault();
            },
            getModifiers: function () {
                return G.modifiers;
            }
        }),
        W = H.extend({
            _class: "KeyEvent",
            initialize: function (t, e, n, i) {
                H.call(this, i),
                    (this.type = t ? "keydown" : "keyup"),
                    (this.key = e),
                    (this.character = n);
            },
            toString: function () {
                return (
                    "{ type: '" +
                    this.type +
                    "', key: '" +
                    this.key +
                    "', character: '" +
                    this.character +
                    "', modifiers: " +
                    this.getModifiers() +
                    " }"
                );
            }
        }),
        G = new (function () {
            function t(t, n, r, h) {
                var u,
                    l = r ? String.fromCharCode(r) : "",
                    c = i[n],
                    d = c || l.toLowerCase(),
                    f = t ? "keydown" : "keyup",
                    _ = Z._focused,
                    g = _ && _.isVisible() && _._scope,
                    p = g && g.tool;
                (o[d] = t),
                    c && (u = e.camelize(c)) in s && (s[u] = t),
                    t ? (a[n] = r) : delete a[n],
                    p &&
                        p.responds(f) &&
                        ((paper = g),
                        p.fire(f, new W(t, d, l, h)),
                        _ && _.update());
            }
            var n,
                i = {
                    8: "backspace",
                    9: "tab",
                    13: "enter",
                    16: "shift",
                    17: "control",
                    18: "option",
                    19: "pause",
                    20: "caps-lock",
                    27: "escape",
                    32: "space",
                    35: "end",
                    36: "home",
                    37: "left",
                    38: "up",
                    39: "right",
                    40: "down",
                    46: "delete",
                    91: "command",
                    93: "command",
                    224: "command"
                },
                r = {
                    9: !0,
                    13: !0,
                    32: !0
                },
                s = new e({
                    shift: !1,
                    control: !1,
                    option: !1,
                    command: !1,
                    capsLock: !1,
                    space: !1
                }),
                a = {},
                o = {};
            return (
                q.add(document, {
                    keydown: function (e) {
                        var a = e.which || e.keyCode;
                        a in i || s.command
                            ? t(!0, a, a in r || s.command ? a : 0, e)
                            : (n = a);
                    },
                    keypress: function (e) {
                        null != n &&
                            (t(!0, n, e.which || e.keyCode, e), (n = null));
                    },
                    keyup: function (e) {
                        var n = e.which || e.keyCode;
                        n in a && t(!1, n, a[n], e);
                    }
                }),
                q.add(window, {
                    blur: function (e) {
                        for (var n in a) t(!1, n, a[n], e);
                    }
                }),
                {
                    modifiers: s,
                    isDown: function (t) {
                        return !!o[t];
                    }
                }
            );
        })(),
        $ = H.extend({
            _class: "MouseEvent",
            initialize: function (t, e, n, i, r) {
                H.call(this, e),
                    (this.type = t),
                    (this.point = n),
                    (this.target = i),
                    (this.delta = r);
            },
            toString: function () {
                return (
                    "{ type: '" +
                    this.type +
                    "', point: " +
                    this.point +
                    ", target: " +
                    this.target +
                    (this.delta ? ", delta: " + this.delta : "") +
                    ", modifiers: " +
                    this.getModifiers() +
                    " }"
                );
            }
        });
    e.extend(n, {
        _class: "Palette",
        _events: ["onChange"],
        initialize: function (n, i, r) {
            var s =
                V.find(".palettejs-panel") ||
                V.find("body").appendChild(
                    V.create("div", {
                        class: "palettejs-panel"
                    })
                );
            (this._element = s.appendChild(
                V.create("table", {
                    class: "palettejs-pane"
                })
            )),
                (this._title = n),
                r || (r = {});
            for (var a in (this.components = i)) {
                var o = i[a];
                o instanceof X ||
                    (null == o.value && (o.value = r[a]),
                    (o.name = a),
                    (o = i[a] = new X(o))),
                    this._element.appendChild(o._element),
                    (o._palette = this),
                    r[a] === t && (r[a] = o.value);
            }
            (this.values = e.each(r, function (t, n) {
                var s = i[n];
                s &&
                    e.define(r, n, {
                        enumerable: !0,
                        configurable: !0,
                        get: function () {
                            return s._value;
                        },
                        set: function (t) {
                            s.setValue(t);
                        }
                    });
            })),
                window.paper && paper.palettes.push(this);
        },
        reset: function () {
            for (var t in this.components) this.components[t].reset();
        },
        remove: function () {
            V.remove(this._element);
        }
    });
    var X = e.extend(n, {
            _class: "Component",
            _events: ["onChange", "onClick"],
            _types: {
                boolean: {
                    type: "checkbox",
                    value: "checked"
                },
                string: {
                    type: "text"
                },
                number: {
                    type: "number",
                    number: !0
                },
                button: {
                    type: "button"
                },
                text: {
                    tag: "div",
                    value: "text"
                },
                slider: {
                    type: "range",
                    number: !0
                },
                list: {
                    tag: "select",
                    setOptions: function () {
                        V.removeChildren(this._input),
                            V.create(
                                e.each(
                                    this._options,
                                    function (t) {
                                        this.push("option", {
                                            value: t,
                                            text: t
                                        });
                                    },
                                    []
                                ),
                                this._input
                            );
                    }
                },
                color: {
                    type: "color",
                    getValue: function (t) {
                        return new D(t);
                    },
                    setValue: function (t) {
                        return new D(t).toCSS(
                            "color" === V.get(this._input, "type")
                        );
                    }
                }
            },
            initialize: function ke(t) {
                (this._id = ke._id = (ke._id || 0) + 1),
                    (this._type =
                        t.type in this._types
                            ? t.type
                            : "options" in t
                            ? "list"
                            : "onClick" in t
                            ? "button"
                            : typeof t.value),
                    (this._meta = this._types[this._type] || {
                        type: this._type
                    });
                var n = this,
                    i = "component-" + this._id;
                (this._dontFire = !0),
                    (this._input = V.create(this._meta.tag || "input", {
                        id: i,
                        type: this._meta.type,
                        events: {
                            change: function () {
                                n.setValue(
                                    V.get(this, n._meta.value || "value")
                                );
                            },
                            click: function () {
                                n.fire("click");
                            }
                        }
                    })),
                    this.attach("change", function (t) {
                        this._dontFire ||
                            this._palette.fire("change", this, this.name, t);
                    }),
                    (this._element = V.create("tr", [
                        "td",
                        [
                            (this._label = V.create("label", {
                                for: i
                            }))
                        ],
                        "td",
                        [this._input]
                    ])),
                    e.each(
                        t,
                        function (t, e) {
                            this[e] = t;
                        },
                        this
                    ),
                    (this._defaultValue = this._value),
                    (this._dontFire = !1);
            },
            getType: function () {
                return this._type;
            },
            getLabel: function () {
                return this.__label;
            },
            setLabel: function (t) {
                (this.__label = t), V.set(this._label, "text", t + ":");
            },
            getOptions: function () {
                return this._options;
            },
            setOptions: function (t) {
                this._options = t;
                var e = this._meta.setOptions;
                e && e.call(this);
            },
            getValue: function () {
                var t = this._value,
                    e = this._meta.getValue;
                return e ? e.call(this, t) : t;
            },
            setValue: function (t) {
                var e = this._meta.value || "value",
                    n = this._meta.setValue;
                n && (t = n.call(this, t)),
                    V.set(this._input, e, t),
                    (t = V.get(this._input, e)),
                    this._meta.number && (t = parseFloat(t, 10)),
                    this._value !== t &&
                        ((this._value = t),
                        this._dontFire || this.fire("change", this.getValue()));
            },
            getRange: function () {
                return [
                    parseFloat(V.get(this._input, "min")),
                    parseFloat(V.get(this._input, "max"))
                ];
            },
            setRange: function (t, e) {
                var n = Array.isArray(t) ? t : [t, e];
                V.set(this._input, {
                    min: n[0],
                    max: n[1]
                });
            },
            getMin: function () {
                return this.getRange()[0];
            },
            setMin: function (t) {
                this.setRange(t, this.getMax());
            },
            getMax: function () {
                return this.getRange()[1];
            },
            setMax: function (t) {
                this.setRange(this.getMin(), t);
            },
            getStep: function () {
                return parseFloat(V.get(this._input, "step"));
            },
            setStep: function (t) {
                V.set(this._input, "step", t);
            },
            reset: function () {
                this.setValue(this._defaultValue);
            }
        }),
        J = H.extend({
            _class: "ToolEvent",
            _item: null,
            initialize: function (t, e, n) {
                (this.tool = t), (this.type = e), (this.event = n);
            },
            _choosePoint: function (t, e) {
                return t ? t : e ? e.clone() : null;
            },
            getPoint: function () {
                return this._choosePoint(this._point, this.tool._point);
            },
            setPoint: function (t) {
                this._point = t;
            },
            getLastPoint: function () {
                return this._choosePoint(this._lastPoint, this.tool._lastPoint);
            },
            setLastPoint: function (t) {
                this._lastPoint = t;
            },
            getDownPoint: function () {
                return this._choosePoint(this._downPoint, this.tool._downPoint);
            },
            setDownPoint: function (t) {
                this._downPoint = t;
            },
            getMiddlePoint: function () {
                return !this._middlePoint && this.tool._lastPoint
                    ? this.tool._point.add(this.tool._lastPoint).divide(2)
                    : this._middlePoint;
            },
            setMiddlePoint: function (t) {
                this._middlePoint = t;
            },
            getDelta: function () {
                return !this._delta && this.tool._lastPoint
                    ? this.tool._point.subtract(this.tool._lastPoint)
                    : this._delta;
            },
            setDelta: function (t) {
                this._delta = t;
            },
            getCount: function () {
                return /^mouse(down|up)$/.test(this.type)
                    ? this.tool._downCount
                    : this.tool._count;
            },
            setCount: function (t) {
                this.tool[
                    /^mouse(down|up)$/.test(this.type) ? "downCount" : "count"
                ] = t;
            },
            getItem: function () {
                if (!this._item) {
                    var t = this.tool._scope.project.hitTest(this.getPoint());
                    if (t) {
                        for (
                            var e = t.item, n = e._parent;
                            /^(Group|CompoundPath)$/.test(n._class);

                        )
                            (e = n), (n = n._parent);
                        this._item = e;
                    }
                }
                return this._item;
            },
            setItem: function (t) {
                this._item = t;
            },
            toString: function () {
                return (
                    "{ type: " +
                    this.type +
                    ", point: " +
                    this.getPoint() +
                    ", count: " +
                    this.getCount() +
                    ", modifiers: " +
                    this.getModifiers() +
                    " }"
                );
            }
        }),
        Y = s.extend({
            _class: "Tool",
            _list: "tools",
            _reference: "tool",
            _events: [
                "onActivate",
                "onDeactivate",
                "onEditOptions",
                "onMouseDown",
                "onMouseUp",
                "onMouseDrag",
                "onMouseMove",
                "onKeyDown",
                "onKeyUp"
            ],
            initialize: function (t) {
                s.call(this),
                    (this._firstMove = !0),
                    (this._count = 0),
                    (this._downCount = 0),
                    this._set(t);
            },
            getMinDistance: function () {
                return this._minDistance;
            },
            setMinDistance: function (t) {
                (this._minDistance = t),
                    null != this._minDistance &&
                        null != this._maxDistance &&
                        this._minDistance > this._maxDistance &&
                        (this._maxDistance = this._minDistance);
            },
            getMaxDistance: function () {
                return this._maxDistance;
            },
            setMaxDistance: function (t) {
                (this._maxDistance = t),
                    null != this._minDistance &&
                        null != this._maxDistance &&
                        this._maxDistance < this._minDistance &&
                        (this._minDistance = t);
            },
            getFixedDistance: function () {
                return this._minDistance == this._maxDistance
                    ? this._minDistance
                    : null;
            },
            setFixedDistance: function (t) {
                (this._minDistance = t), (this._maxDistance = t);
            },
            _updateEvent: function (t, e, n, i, r, s, a) {
                if (!r) {
                    if (null != n || null != i) {
                        var o = null != n ? n : 0,
                            h = e.subtract(this._point),
                            u = h.getLength();
                        if (o > u) return !1;
                        var l = null != i ? i : 0;
                        if (0 != l)
                            if (u > l) e = this._point.add(h.normalize(l));
                            else if (a) return !1;
                    }
                    if (s && e.equals(this._point)) return !1;
                }
                switch (
                    ((this._lastPoint =
                        r && "mousemove" == t ? e : this._point),
                    (this._point = e),
                    t)
                ) {
                    case "mousedown":
                        (this._lastPoint = this._downPoint),
                            (this._downPoint = this._point),
                            this._downCount++;
                        break;
                    case "mouseup":
                        this._lastPoint = this._downPoint;
                }
                return (this._count = r ? 0 : this._count + 1), !0;
            },
            _fireEvent: function (t, e) {
                var n = paper.project._removeSets;
                if (n) {
                    "mouseup" === t && (n.mousedrag = null);
                    var i = n[t];
                    if (i) {
                        for (var r in i) {
                            var s = i[r];
                            for (var a in n) {
                                var o = n[a];
                                o && o != i && delete o[s._id];
                            }
                            s.remove();
                        }
                        n[t] = null;
                    }
                }
                return this.responds(t) && this.fire(t, new J(this, t, e));
            },
            _handleEvent: function (t, e, n) {
                paper = this._scope;
                var i = !1;
                switch (t) {
                    case "mousedown":
                        this._updateEvent(t, e, null, null, !0, !1, !1),
                            (i = this._fireEvent(t, n));
                        break;
                    case "mousedrag":
                        for (
                            var r = !1, s = !1;
                            this._updateEvent(
                                t,
                                e,
                                this.minDistance,
                                this.maxDistance,
                                !1,
                                r,
                                s
                            );

                        )
                            (i = this._fireEvent(t, n) || i),
                                (r = !0),
                                (s = !0);
                        break;
                    case "mouseup":
                        !e.equals(this._point) &&
                            this._updateEvent(
                                "mousedrag",
                                e,
                                this.minDistance,
                                this.maxDistance,
                                !1,
                                !1,
                                !1
                            ) &&
                            (i = this._fireEvent("mousedrag", n)),
                            this._updateEvent(
                                t,
                                e,
                                null,
                                this.maxDistance,
                                !1,
                                !1,
                                !1
                            ),
                            (i = this._fireEvent(t, n) || i),
                            this._updateEvent(t, e, null, null, !0, !1, !1),
                            (this._firstMove = !0);
                        break;
                    case "mousemove":
                        for (
                            ;
                            this._updateEvent(
                                t,
                                e,
                                this.minDistance,
                                this.maxDistance,
                                this._firstMove,
                                !0,
                                !1
                            );

                        )
                            (i = this._fireEvent(t, n) || i),
                                (this._firstMove = !1);
                }
                return i && n.preventDefault(), i;
            }
        }),
        K = {
            request: function (t, e, n) {
                var i = new (window.ActiveXObject || XMLHttpRequest)(
                    "Microsoft.XMLHTTP"
                );
                return (
                    i.open(t.toUpperCase(), e, !0),
                    "overrideMimeType" in i && i.overrideMimeType("text/plain"),
                    (i.onreadystatechange = function () {
                        if (4 === i.readyState) {
                            var t = i.status;
                            if (0 !== t && 200 !== t)
                                throw Error(
                                    "Could not load " + e + " (Error " + t + ")"
                                );
                            n.call(i, i.responseText);
                        }
                    }),
                    i.send(null)
                );
            }
        },
        Q = {
            canvases: [],
            getCanvas: function (t, e) {
                var n,
                    i = !0;
                "object" == typeof t && ((e = t.height), (t = t.width)),
                    (n = this.canvases.length
                        ? this.canvases.pop()
                        : document.createElement("canvas"));
                var r = n.getContext("2d");
                return (
                    n.width === t && n.height === e
                        ? i && r.clearRect(0, 0, t + 1, e + 1)
                        : ((n.width = t), (n.height = e)),
                    r.save(),
                    n
                );
            },
            getContext: function (t, e) {
                return this.getCanvas(t, e).getContext("2d");
            },
            release: function (t) {
                var e = t.canvas ? t.canvas : t;
                e.getContext("2d").restore(), this.canvases.push(e);
            }
        },
        te = new (function () {
            function t(t, e, n) {
                return 0.2989 * t + 0.587 * e + 0.114 * n;
            }
            function n(e, n, i, r) {
                var s = r - t(e, n, i);
                (f = e + s), (_ = n + s), (g = i + s);
                var r = t(f, _, g),
                    a = p(f, _, g),
                    o = v(f, _, g);
                if (0 > a) {
                    var h = r - a;
                    (f = r + ((f - r) * r) / h),
                        (_ = r + ((_ - r) * r) / h),
                        (g = r + ((g - r) * r) / h);
                }
                if (o > 255) {
                    var u = 255 - r,
                        l = o - r;
                    (f = r + ((f - r) * u) / l),
                        (_ = r + ((_ - r) * u) / l),
                        (g = r + ((g - r) * u) / l);
                }
            }
            function i(t, e, n) {
                return v(t, e, n) - p(t, e, n);
            }
            function r(t, e, n, i) {
                var r,
                    s = [t, e, n],
                    a = v(t, e, n),
                    o = p(t, e, n);
                (o = o === t ? 0 : o === e ? 1 : 2),
                    (a = a === t ? 0 : a === e ? 1 : 2),
                    (r = 0 === p(o, a) ? (1 === v(o, a) ? 2 : 1) : 0),
                    s[a] > s[o]
                        ? ((s[r] = ((s[r] - s[o]) * i) / (s[a] - s[o])),
                          (s[a] = i))
                        : (s[r] = s[a] = 0),
                    (s[o] = 0),
                    (f = s[0]),
                    (_ = s[1]),
                    (g = s[2]);
            }
            var s,
                a,
                o,
                h,
                u,
                l,
                c,
                d,
                f,
                _,
                g,
                p = Math.min,
                v = Math.max,
                m = Math.abs,
                y = {
                    multiply: function () {
                        (f = (u * s) / 255),
                            (_ = (l * a) / 255),
                            (g = (c * o) / 255);
                    },
                    screen: function () {
                        (f = u + s - (u * s) / 255),
                            (_ = l + a - (l * a) / 255),
                            (g = c + o - (c * o) / 255);
                    },
                    overlay: function () {
                        (f =
                            128 > u
                                ? (2 * u * s) / 255
                                : 255 - (2 * (255 - u) * (255 - s)) / 255),
                            (_ =
                                128 > l
                                    ? (2 * l * a) / 255
                                    : 255 - (2 * (255 - l) * (255 - a)) / 255),
                            (g =
                                128 > c
                                    ? (2 * c * o) / 255
                                    : 255 - (2 * (255 - c) * (255 - o)) / 255);
                    },
                    "soft-light": function () {
                        var t = (s * u) / 255;
                        (f =
                            t +
                            (u * (255 - ((255 - u) * (255 - s)) / 255 - t)) /
                                255),
                            (t = (a * l) / 255),
                            (_ =
                                t +
                                (l *
                                    (255 - ((255 - l) * (255 - a)) / 255 - t)) /
                                    255),
                            (t = (o * c) / 255),
                            (g =
                                t +
                                (c *
                                    (255 - ((255 - c) * (255 - o)) / 255 - t)) /
                                    255);
                    },
                    "hard-light": function () {
                        (f =
                            128 > s
                                ? (2 * s * u) / 255
                                : 255 - (2 * (255 - s) * (255 - u)) / 255),
                            (_ =
                                128 > a
                                    ? (2 * a * l) / 255
                                    : 255 - (2 * (255 - a) * (255 - l)) / 255),
                            (g =
                                128 > o
                                    ? (2 * o * c) / 255
                                    : 255 - (2 * (255 - o) * (255 - c)) / 255);
                    },
                    "color-dodge": function () {
                        (f =
                            0 === u
                                ? 0
                                : 255 === s
                                ? 255
                                : p(255, (255 * u) / (255 - s))),
                            (_ =
                                0 === l
                                    ? 0
                                    : 255 === a
                                    ? 255
                                    : p(255, (255 * l) / (255 - a))),
                            (g =
                                0 === c
                                    ? 0
                                    : 255 === o
                                    ? 255
                                    : p(255, (255 * c) / (255 - o)));
                    },
                    "color-burn": function () {
                        (f =
                            255 === u
                                ? 255
                                : 0 === s
                                ? 0
                                : v(0, 255 - (255 * (255 - u)) / s)),
                            (_ =
                                255 === l
                                    ? 255
                                    : 0 === a
                                    ? 0
                                    : v(0, 255 - (255 * (255 - l)) / a)),
                            (g =
                                255 === c
                                    ? 255
                                    : 0 === o
                                    ? 0
                                    : v(0, 255 - (255 * (255 - c)) / o));
                    },
                    darken: function () {
                        (f = s > u ? u : s),
                            (_ = a > l ? l : a),
                            (g = o > c ? c : o);
                    },
                    lighten: function () {
                        (f = u > s ? u : s),
                            (_ = l > a ? l : a),
                            (g = c > o ? c : o);
                    },
                    difference: function () {
                        (f = u - s),
                            0 > f && (f = -f),
                            (_ = l - a),
                            0 > _ && (_ = -_),
                            (g = c - o),
                            0 > g && (g = -g);
                    },
                    exclusion: function () {
                        (f = u + (s * (255 - u - u)) / 255),
                            (_ = l + (a * (255 - l - l)) / 255),
                            (g = c + (o * (255 - c - c)) / 255);
                    },
                    hue: function () {
                        r(s, a, o, i(u, l, c)), n(f, _, g, t(u, l, c));
                    },
                    saturation: function () {
                        r(u, l, c, i(s, a, o)), n(f, _, g, t(u, l, c));
                    },
                    luminosity: function () {
                        n(u, l, c, t(s, a, o));
                    },
                    color: function () {
                        n(s, a, o, t(u, l, c));
                    },
                    add: function () {
                        (f = p(u + s, 255)),
                            (_ = p(l + a, 255)),
                            (g = p(c + o, 255));
                    },
                    subtract: function () {
                        (f = v(u - s, 0)), (_ = v(l - a, 0)), (g = v(c - o, 0));
                    },
                    average: function () {
                        (f = (u + s) / 2), (_ = (l + a) / 2), (g = (c + o) / 2);
                    },
                    negation: function () {
                        (f = 255 - m(255 - s - u)),
                            (_ = 255 - m(255 - a - l)),
                            (g = 255 - m(255 - o - c));
                    }
                },
                w = (this.nativeModes = e.each(
                    [
                        "source-over",
                        "source-in",
                        "source-out",
                        "source-atop",
                        "destination-over",
                        "destination-in",
                        "destination-out",
                        "destination-atop",
                        "lighter",
                        "darker",
                        "copy",
                        "xor"
                    ],
                    function (t) {
                        this[t] = !0;
                    },
                    {}
                )),
                x = Q.getContext(1, 1);
            e.each(y, function (t, e) {
                var n = "darken" === e,
                    i = !1;
                x.save();
                try {
                    (x.fillStyle = n ? "#300" : "#a00"),
                        x.fillRect(0, 0, 1, 1),
                        (x.globalCompositeOperation = e),
                        x.globalCompositeOperation === e &&
                            ((x.fillStyle = n ? "#a00" : "#300"),
                            x.fillRect(0, 0, 1, 1),
                            (i =
                                x.getImageData(0, 0, 1, 1).data[0] !== n
                                    ? 170
                                    : 51));
                } catch (r) {}
                x.restore(), (w[e] = i);
            }),
                Q.release(x),
                (this.process = function (t, e, n, i, r) {
                    var p = e.canvas,
                        v = "normal" === t;
                    if (v || w[t])
                        n.save(),
                            n.setTransform(1, 0, 0, 1, 0, 0),
                            (n.globalAlpha = i),
                            v || (n.globalCompositeOperation = t),
                            n.drawImage(p, r.x, r.y),
                            n.restore();
                    else {
                        var m = y[t];
                        if (!m) return;
                        for (
                            var x = n.getImageData(r.x, r.y, p.width, p.height),
                                b = x.data,
                                C = e.getImageData(
                                    0,
                                    0,
                                    p.width,
                                    p.height
                                ).data,
                                S = 0,
                                P = b.length;
                            P > S;
                            S += 4
                        ) {
                            (s = C[S]),
                                (u = b[S]),
                                (a = C[S + 1]),
                                (l = b[S + 1]),
                                (o = C[S + 2]),
                                (c = b[S + 2]),
                                (h = C[S + 3]),
                                (d = b[S + 3]),
                                m();
                            var k = (h * i) / 255,
                                M = 1 - k;
                            (b[S] = k * f + M * u),
                                (b[S + 1] = k * _ + M * l),
                                (b[S + 2] = k * g + M * c),
                                (b[S + 3] = h * i + M * d);
                        }
                        n.putImageData(x, r.x, r.y);
                    }
                });
        })(),
        ee = e.each(
            {
                fillColor: ["fill", "color"],
                strokeColor: ["stroke", "color"],
                strokeWidth: ["stroke-width", "number"],
                strokeCap: ["stroke-linecap", "string"],
                strokeJoin: ["stroke-linejoin", "string"],
                strokeScaling: [
                    "vector-effect",
                    "lookup",
                    {
                        true: "none",
                        false: "non-scaling-stroke"
                    },
                    function (t, e) {
                        return (
                            !e &&
                            (t instanceof I || t instanceof b || t instanceof N)
                        );
                    }
                ],
                miterLimit: ["stroke-miterlimit", "number"],
                dashArray: ["stroke-dasharray", "array"],
                dashOffset: ["stroke-dashoffset", "number"],
                fontFamily: ["font-family", "string"],
                fontWeight: ["font-weight", "string"],
                fontSize: ["font-size", "number"],
                justification: [
                    "text-anchor",
                    "lookup",
                    {
                        left: "start",
                        center: "middle",
                        right: "end"
                    }
                ],
                opacity: ["opacity", "number"],
                blendMode: ["mix-blend-mode", "string"]
            },
            function (t, n) {
                var i = e.capitalize(n),
                    r = t[2];
                this[n] = {
                    type: t[1],
                    property: n,
                    attribute: t[0],
                    toSVG: r,
                    fromSVG:
                        r &&
                        e.each(
                            r,
                            function (t, e) {
                                this[t] = e;
                            },
                            {}
                        ),
                    exportFilter: t[3],
                    get: "get" + i,
                    set: "set" + i
                };
            },
            {}
        ),
        ne = {
            href: "http://www.w3.org/1999/xlink",
            xlink: "http://www.w3.org/2000/xmlns"
        };
    return (
        new (function () {
            function t(t, e) {
                for (var n in e) {
                    var i = e[n],
                        r = ne[n];
                    "number" == typeof i && (i = S.number(i)),
                        r ? t.setAttributeNS(r, n, i) : t.setAttribute(n, i);
                }
                return t;
            }
            function n(e, n) {
                return t(
                    document.createElementNS("http://www.w3.org/2000/svg", e),
                    n
                );
            }
            function r(t, n, i) {
                var r = new e(),
                    s = t.getTranslation();
                if (n) {
                    t = t.shiftless();
                    var a = t._inverseTransform(s);
                    (r[i ? "cx" : "x"] = a.x),
                        (r[i ? "cy" : "y"] = a.y),
                        (s = null);
                }
                if (!t.isIdentity()) {
                    var h = t.decompose();
                    if (h && !h.shearing) {
                        var u = [],
                            l = h.rotation,
                            c = h.scaling;
                        s &&
                            !s.isZero() &&
                            u.push("translate(" + S.point(s) + ")"),
                            l && u.push("rotate(" + S.number(l) + ")"),
                            (o.isZero(c.x - 1) && o.isZero(c.y - 1)) ||
                                u.push("scale(" + S.point(c) + ")"),
                            (r.transform = u.join(" "));
                    } else
                        r.transform = "matrix(" + t.getValues().join(",") + ")";
                }
                return r;
            }
            function s(e, i) {
                for (
                    var s = r(e._matrix),
                        a = e._children,
                        o = n("g", s),
                        h = 0,
                        u = a.length;
                    u > h;
                    h++
                ) {
                    var l = a[h],
                        c = b(l, i);
                    if (c)
                        if (l.isClipMask()) {
                            var d = n("clipPath");
                            d.appendChild(c),
                                w(l, d, "clip"),
                                t(o, {
                                    "clip-path": "url(#" + d.id + ")"
                                });
                        } else o.appendChild(c);
                }
                return o;
            }
            function h(t) {
                var e = r(t._matrix, !0),
                    i = t.getSize();
                return (
                    (e.x -= i.width / 2),
                    (e.y -= i.height / 2),
                    (e.width = i.width),
                    (e.height = i.height),
                    (e.href = t.toDataURL()),
                    n("image", e)
                );
            }
            function u(t, e) {
                if (e.matchShapes) {
                    var s = t.toShape(!1);
                    if (s) return c(s, e);
                }
                var a,
                    o = t._segments,
                    h = r(t._matrix);
                if (0 === o.length) return null;
                if (t.isPolygon())
                    if (o.length >= 3) {
                        a = t._closed ? "polygon" : "polyline";
                        var u = [];
                        for (i = 0, l = o.length; l > i; i++)
                            u.push(S.point(o[i]._point));
                        h.points = u.join(" ");
                    } else {
                        a = "line";
                        var d = o[0]._point,
                            f = o[o.length - 1]._point;
                        h.set({
                            x1: d.x,
                            y1: d.y,
                            x2: f.x,
                            y2: f.y
                        });
                    }
                else (a = "path"), (h.d = t.getPathData(null, e.precision));
                return n(a, h);
            }
            function c(t) {
                var e = t._type,
                    i = t._radius,
                    s = r(t._matrix, !0, "rectangle" !== e);
                if ("rectangle" === e) {
                    e = "rect";
                    var a = t._size,
                        o = a.width,
                        h = a.height;
                    (s.x -= o / 2),
                        (s.y -= h / 2),
                        (s.width = o),
                        (s.height = h),
                        i.isZero() && (i = null);
                }
                return (
                    i &&
                        ("circle" === e
                            ? (s.r = i)
                            : ((s.rx = i.width), (s.ry = i.height))),
                    n(e, s)
                );
            }
            function d(t, e) {
                var i = r(t._matrix),
                    s = t.getPathData(null, e.precision);
                return s && (i.d = s), n("path", i);
            }
            function f(t, e) {
                var i = r(t._matrix, !0),
                    s = t.getSymbol(),
                    a = m(s, "symbol"),
                    o = s.getDefinition(),
                    h = o.getBounds();
                return (
                    a ||
                        ((a = n("symbol", {
                            viewBox: S.rectangle(h)
                        })),
                        a.appendChild(b(o, e)),
                        w(s, a, "symbol")),
                    (i.href = "#" + a.id),
                    (i.x += h.x),
                    (i.y += h.y),
                    (i.width = S.number(h.width)),
                    (i.height = S.number(h.height)),
                    n("use", i)
                );
            }
            function _(t) {
                var e = m(t, "color");
                if (!e) {
                    var i,
                        r = t.getGradient(),
                        s = r._radial,
                        a = t.getOrigin().transform(),
                        o = t.getDestination().transform();
                    if (s) {
                        i = {
                            cx: a.x,
                            cy: a.y,
                            r: a.getDistance(o)
                        };
                        var h = t.getHighlight();
                        h && ((h = h.transform()), (i.fx = h.x), (i.fy = h.y));
                    } else
                        i = {
                            x1: a.x,
                            y1: a.y,
                            x2: o.x,
                            y2: o.y
                        };
                    (i.gradientUnits = "userSpaceOnUse"),
                        (e = n((s ? "radial" : "linear") + "Gradient", i));
                    for (var u = r._stops, l = 0, c = u.length; c > l; l++) {
                        var d = u[l],
                            f = d._color,
                            _ = f.getAlpha();
                        (i = {
                            offset: d._rampPoint,
                            "stop-color": f.toCSS(!0)
                        }),
                            1 > _ && (i["stop-opacity"] = _),
                            e.appendChild(n("stop", i));
                    }
                    w(t, e, "color");
                }
                return "url(#" + e.id + ")";
            }
            function g(t) {
                var e = n("text", r(t._matrix, !0));
                return (e.textContent = t._content), e;
            }
            function p(n, i, r) {
                var s = {},
                    a = !r && n.getParent();
                return (
                    null != n._name && (s.id = n._name),
                    e.each(ee, function (t) {
                        var i = t.get,
                            r = t.type,
                            o = n[i]();
                        if (
                            t.exportFilter
                                ? t.exportFilter(n, o)
                                : !a || !e.equals(a[i](), o)
                        ) {
                            if ("color" === r && null != o) {
                                var h = o.getAlpha();
                                1 > h && (s[t.attribute + "-opacity"] = h);
                            }
                            s[t.attribute] =
                                null == o
                                    ? "none"
                                    : "number" === r
                                    ? S.number(o)
                                    : "color" === r
                                    ? o.gradient
                                        ? _(o, n)
                                        : o.toCSS(!0)
                                    : "array" === r
                                    ? o.join(",")
                                    : "lookup" === r
                                    ? t.toSVG[o]
                                    : o;
                        }
                    }),
                    1 === s.opacity && delete s.opacity,
                    n._visible || (s.visibility = "hidden"),
                    t(i, s)
                );
            }
            function m(t, e) {
                return (
                    P ||
                        (P = {
                            ids: {},
                            svgs: {}
                        }),
                    t && P.svgs[e + "-" + t._id]
                );
            }
            function w(t, e, n) {
                P || m();
                var i = (P.ids[n] = (P.ids[n] || 0) + 1);
                (e.id = n + "-" + i), (P.svgs[n + "-" + t._id] = e);
            }
            function x(t, e) {
                var i = t,
                    r = null;
                if (P) {
                    i = "svg" === t.nodeName.toLowerCase() && t;
                    for (var s in P.svgs)
                        r ||
                            (i || ((i = n("svg")), i.appendChild(t)),
                            (r = i.insertBefore(n("defs"), i.firstChild))),
                            r.appendChild(P.svgs[s]);
                    P = null;
                }
                return e.asString
                    ? new XMLSerializer().serializeToString(i)
                    : i;
            }
            function b(t, e, n) {
                var i = k[t._class],
                    r = i && i(t, e);
                if (r) {
                    var s = e.onExport;
                    s && (r = s(t, r, e) || r);
                    var a = JSON.stringify(t._data);
                    a && "{}" !== a && r.setAttribute("data-paper-data", a);
                }
                return r && p(t, r, n);
            }
            function C(t) {
                return t || (t = {}), (S = new a(t.precision)), t;
            }
            var S,
                P,
                k = {
                    Group: s,
                    Layer: s,
                    Raster: h,
                    Path: u,
                    Shape: c,
                    CompoundPath: d,
                    PlacedSymbol: f,
                    PointText: g
                };
            y.inject({
                exportSVG: function (t) {
                    return (t = C(t)), x(b(this, t, !0), t);
                }
            }),
                v.inject({
                    exportSVG: function (t) {
                        t = C(t);
                        var e = this.layers,
                            i = this.getView(),
                            s = i.getViewSize(),
                            a = n("svg", {
                                x: 0,
                                y: 0,
                                width: s.width,
                                height: s.height,
                                version: "1.1",
                                xmlns: "http://www.w3.org/2000/svg",
                                "xmlns:xlink": "http://www.w3.org/1999/xlink"
                            }),
                            o = a,
                            h = i._matrix;
                        h.isIdentity() || (o = a.appendChild(n("g", r(h))));
                        for (var u = 0, l = e.length; l > u; u++)
                            o.appendChild(b(e[u], t, !0));
                        return x(a, t);
                    }
                });
        })(),
        new (function () {
            function n(t, e, n, i) {
                var r = ne[e],
                    s = r ? t.getAttributeNS(r, e) : t.getAttribute(e);
                return (
                    "null" === s && (s = null),
                    null == s ? (i ? null : n ? "" : 0) : n ? s : parseFloat(s)
                );
            }
            function i(t, e, i, r) {
                return (
                    (e = n(t, e, !1, r)),
                    (i = n(t, i, !1, r)),
                    !r || (null != e && null != i) ? new h(e, i) : null
                );
            }
            function r(t, e, i, r) {
                return (
                    (e = n(t, e, !1, r)),
                    (i = n(t, i, !1, r)),
                    !r || (null != e && null != i) ? new c(e, i) : null
                );
            }
            function s(t, e, n) {
                return "none" === t
                    ? null
                    : "number" === e
                    ? parseFloat(t)
                    : "array" === e
                    ? t
                        ? t.split(/[\s,]+/g).map(parseFloat)
                        : []
                    : "color" === e
                    ? S(t) || t
                    : "lookup" === e
                    ? n[t]
                    : t;
            }
            function a(t, e, n, i) {
                var r = t.childNodes,
                    s = "clippath" === e,
                    a = new w(),
                    o = a._project,
                    h = o._currentStyle,
                    u = [];
                s || ((a = x(a, t, i)), (o._currentStyle = a._style.clone()));
                for (var l = 0, c = r.length; c > l; l++) {
                    var d,
                        f = r[l];
                    1 !== f.nodeType ||
                        !(d = P(f, n, !1)) ||
                        d instanceof m ||
                        u.push(d);
                }
                return (
                    a.addChildren(u),
                    s && (a = x(a.reduce(), t, i)),
                    (o._currentStyle = h),
                    (s || "defs" === e) && (a.remove(), (a = null)),
                    a
                );
            }
            function o(t, e) {
                for (
                    var n = t
                            .getAttribute("points")
                            .match(
                                /[+-]?(?:\d*\.\d+|\d+\.?)(?:[eE][+-]?\d+)?/g
                            ),
                        i = [],
                        r = 0,
                        s = n.length;
                    s > r;
                    r += 2
                )
                    i.push(new h(parseFloat(n[r]), parseFloat(n[r + 1])));
                var a = new O(i);
                return "polygon" === e && a.closePath(), a;
            }
            function u(t) {
                var e = t.getAttribute("d"),
                    n = {
                        pathData: e
                    };
                return e.match(/m/gi).length > 1 || /z\S+/i.test(e)
                    ? new T(n)
                    : new O(n);
            }
            function l(t, e) {
                var r,
                    s = (n(t, "href", !0) || "").substring(1),
                    a = "radialgradient" === e;
                if (s) r = z[s].getGradient();
                else {
                    for (
                        var o = t.childNodes, h = [], u = 0, l = o.length;
                        l > u;
                        u++
                    ) {
                        var c = o[u];
                        1 === c.nodeType && h.push(x(new R(), c));
                    }
                    r = new B(h, a);
                }
                var d, f, _;
                return (
                    a
                        ? ((d = i(t, "cx", "cy")),
                          (f = d.add(n(t, "r"), 0)),
                          (_ = i(t, "fx", "fy", !0)))
                        : ((d = i(t, "x1", "y1")), (f = i(t, "x2", "y2"))),
                    x(new D(r, d, f, _), t),
                    null
                );
            }
            function d(t, e, n, i) {
                for (
                    var r = (i.getAttribute(n) || "").split(/\)\s*/g),
                        s = new g(),
                        a = 0,
                        o = r.length;
                    o > a;
                    a++
                ) {
                    var h = r[a];
                    if (!h) break;
                    for (
                        var u = h.split("("),
                            l = u[0],
                            c = u[1].split(/[\s,]+/g),
                            d = 0,
                            f = c.length;
                        f > d;
                        d++
                    )
                        c[d] = parseFloat(c[d]);
                    switch (l) {
                        case "matrix":
                            s.concatenate(
                                new g(c[0], c[1], c[2], c[3], c[4], c[5])
                            );
                            break;
                        case "rotate":
                            s.rotate(c[0], c[1], c[2]);
                            break;
                        case "translate":
                            s.translate(c[0], c[1]);
                            break;
                        case "scale":
                            s.scale(c);
                            break;
                        case "skewX":
                            s.skew(c[0], 0);
                            break;
                        case "skewY":
                            s.skew(0, c[0]);
                    }
                }
                t.transform(s);
            }
            function _(t, e, n) {
                var i =
                    t[
                        "fill-opacity" === n ? "getFillColor" : "getStrokeColor"
                    ]();
                i && i.setAlpha(parseFloat(e));
            }
            function p(n, i, r) {
                var s = n.attributes[i],
                    a = s && s.value;
                if (!a) {
                    var o = e.camelize(i);
                    (a = n.style[o]),
                        a || r.node[o] === r.parent[o] || (a = r.node[o]);
                }
                return a ? ("none" === a ? null : a) : t;
            }
            function x(n, i, r) {
                var s = {
                    node: V.getStyles(i) || {},
                    parent: (!r && V.getStyles(i.parentNode)) || {}
                };
                return (
                    e.each(M, function (r, a) {
                        var o = p(i, a, s);
                        o !== t && (n = e.pick(r(n, o, a, i, s), n));
                    }),
                    n
                );
            }
            function S(t) {
                var e = t && t.match(/\((?:#|)([^)']+)/);
                return e && z[e[1]];
            }
            function P(t, e, n) {
                function i(t) {
                    paper = s;
                    var i = P(t, e, n),
                        r = e.onLoad,
                        a = s.project && s.getView();
                    r && r.call(this, i), a.update();
                }
                if (!t) return null;
                e
                    ? "function" == typeof e &&
                      (e = {
                          onLoad: e
                      })
                    : (e = {});
                var r = t,
                    s = paper;
                if (n)
                    if ("string" != typeof t || /^.*</.test(t)) {
                        if ("undefined" != typeof File && t instanceof File) {
                            var a = new FileReader();
                            return (
                                (a.onload = function () {
                                    i(a.result);
                                }),
                                a.readAsText(t)
                            );
                        }
                    } else {
                        if (((r = document.getElementById(t)), !r))
                            return K.request("get", t, i);
                        t = null;
                    }
                if (
                    ("string" == typeof t &&
                        (r = new DOMParser().parseFromString(
                            t,
                            "image/svg+xml"
                        )),
                    !r.nodeName)
                )
                    throw Error("Unsupported SVG source: " + t);
                var o,
                    h = r.nodeName.toLowerCase(),
                    u = k[h],
                    l = r.getAttribute && r.getAttribute("data-paper-data"),
                    c = s.settings,
                    d = c.applyMatrix;
                if (
                    ((c.applyMatrix = !1),
                    (o = (u && u(r, h, e, n)) || null),
                    (c.applyMatrix = d),
                    o)
                ) {
                    "#document" === h || o instanceof w || (o = x(o, r, n));
                    var f = e.onImport;
                    f && (o = f(r, o, e) || o),
                        e.expandShapes &&
                            o instanceof b &&
                            (o.remove(), (o = o.toPath())),
                        l && (o._data = JSON.parse(l));
                }
                return n && (z = {}), o;
            }
            var k = {
                    "#document": function (t, e, n, i) {
                        for (
                            var r = t.childNodes, s = 0, a = r.length;
                            a > s;
                            s++
                        ) {
                            var o = r[s];
                            if (1 === o.nodeType) {
                                var h = o.nextSibling;
                                document.body.appendChild(o);
                                var u = P(o, n, i);
                                return (
                                    h ? t.insertBefore(o, h) : t.appendChild(o),
                                    u
                                );
                            }
                        }
                    },
                    g: a,
                    svg: a,
                    clippath: a,
                    polygon: o,
                    polyline: o,
                    path: u,
                    lineargradient: l,
                    radialgradient: l,
                    image: function (t) {
                        var e = new C(n(t, "href", !0));
                        return (
                            e.attach("load", function () {
                                var e = r(t, "width", "height");
                                this.setSize(e);
                                var n = this._matrix._transformPoint(
                                    i(t, "x", "y").add(e.divide(2))
                                );
                                this.translate(n);
                            }),
                            e
                        );
                    },
                    symbol: function (t, e, n, i) {
                        return new m(a(t, e, n, i), !0);
                    },
                    defs: a,
                    use: function (t) {
                        var e = (n(t, "href", !0) || "").substring(1),
                            r = z[e],
                            s = i(t, "x", "y");
                        return r
                            ? r instanceof m
                                ? r.place(s)
                                : r.clone().translate(s)
                            : null;
                    },
                    circle: function (t) {
                        return new b.Circle(i(t, "cx", "cy"), n(t, "r"));
                    },
                    ellipse: function (t) {
                        return new b.Ellipse({
                            center: i(t, "cx", "cy"),
                            radius: r(t, "rx", "ry")
                        });
                    },
                    rect: function (t) {
                        var e = i(t, "x", "y"),
                            n = r(t, "width", "height"),
                            s = r(t, "rx", "ry");
                        return new b.Rectangle(new f(e, n), s);
                    },
                    line: function (t) {
                        return new O.Line(i(t, "x1", "y1"), i(t, "x2", "y2"));
                    },
                    text: function (t) {
                        var e = new j(i(t, "x", "y").add(i(t, "dx", "dy")));
                        return e.setContent(t.textContent.trim() || ""), e;
                    }
                },
                M = e.each(
                    ee,
                    function (t) {
                        this[t.attribute] = function (e, n) {
                            if (
                                (e[t.set](s(n, t.type, t.fromSVG)),
                                "color" === t.type && e instanceof b)
                            ) {
                                var i = e[t.get]();
                                i &&
                                    i.transform(
                                        new g().translate(
                                            e.getPosition(!0).negate()
                                        )
                                    );
                            }
                        };
                    },
                    {
                        id: function (t, e) {
                            (z[e] = t), t.setName && t.setName(e);
                        },
                        "clip-path": function (t, e) {
                            var n = S(e);
                            if (n) {
                                if (
                                    ((n = n.clone()),
                                    n.setClipMask(!0),
                                    !(t instanceof w))
                                )
                                    return new w(n, t);
                                t.insertChild(0, n);
                            }
                        },
                        gradientTransform: d,
                        transform: d,
                        "fill-opacity": _,
                        "stroke-opacity": _,
                        visibility: function (t, e) {
                            t.setVisible("visible" === e);
                        },
                        display: function (t, e) {
                            t.setVisible(null !== e);
                        },
                        "stop-color": function (t, e) {
                            t.setColor && t.setColor(e);
                        },
                        "stop-opacity": function (t, e) {
                            t._color && t._color.setAlpha(parseFloat(e));
                        },
                        offset: function (t, e) {
                            var n = e.match(/(.*)%$/);
                            t.setRampPoint(n ? n[1] / 100 : parseFloat(e));
                        },
                        viewBox: function (t, e, n, i, a) {
                            var o = new f(s(e, "array")),
                                h = r(i, "width", "height", !0);
                            if (t instanceof w) {
                                var u = h ? o.getSize().divide(h) : 1,
                                    l = new g()
                                        .translate(o.getPoint())
                                        .scale(u);
                                t.transform(l.inverted());
                            } else if (t instanceof m) {
                                h && o.setSize(h);
                                var c = "visible" != p(i, "overflow", a),
                                    d = t._definition;
                                c &&
                                    !o.contains(d.getBounds()) &&
                                    ((c = new b.Rectangle(o).transform(
                                        d._matrix
                                    )),
                                    c.setClipMask(!0),
                                    d.addChild(c));
                            }
                        }
                    }
                ),
                z = {};
            y.inject({
                importSVG: function (t, e) {
                    return this.addChild(P(t, e, !0));
                }
            }),
                v.inject({
                    importSVG: function (t, e) {
                        return this.activate(), P(t, e, !0);
                    }
                });
        })(),
        (e.exports.PaperScript = function () {
            function t(t, e, n) {
                var i = v[e];
                if (t && t[i]) {
                    var r = t[i](n);
                    return "!=" === e ? !r : r;
                }
                switch (e) {
                    case "+":
                        return t + n;
                    case "-":
                        return t - n;
                    case "*":
                        return t * n;
                    case "/":
                        return t / n;
                    case "%":
                        return t % n;
                    case "==":
                        return t == n;
                    case "!=":
                        return t != n;
                }
            }
            function n(t, e) {
                var n = m[t];
                if (n && e && e[n]) return e[n]();
                switch (t) {
                    case "+":
                        return +e;
                    case "-":
                        return -e;
                }
            }
            function i(t, e) {
                return _.acorn.parse(t, e);
            }
            function s(t, e, n) {
                function r(t) {
                    for (var e = 0, n = h.length; n > e; e++) {
                        var i = h[e];
                        if (i[0] >= t) break;
                        t += i[1];
                    }
                    return t;
                }
                function s(e) {
                    return t.substring(r(e.range[0]), r(e.range[1]));
                }
                function a(e, n) {
                    for (
                        var i = r(e.range[0]),
                            s = r(e.range[1]),
                            a = 0,
                            o = h.length - 1;
                        o >= 0;
                        o--
                    )
                        if (i > h[o][0]) {
                            a = o + 1;
                            break;
                        }
                    h.splice(a, 0, [i, n.length - s + i]),
                        (t = t.substring(0, i) + n + t.substring(s));
                }
                function o(t, e) {
                    if (t) {
                        for (var n in t)
                            if ("range" !== n && "loc" !== n) {
                                var i = t[n];
                                if (Array.isArray(i))
                                    for (var r = 0, h = i.length; h > r; r++)
                                        o(i[r], t);
                                else i && "object" == typeof i && o(i, t);
                            }
                        switch (t.type) {
                            case "UnaryExpression":
                                if (
                                    t.operator in m &&
                                    "Literal" !== t.argument.type
                                ) {
                                    var u = s(t.argument);
                                    a(
                                        t,
                                        '$__("' + t.operator + '", ' + u + ")"
                                    );
                                }
                                break;
                            case "BinaryExpression":
                                if (
                                    t.operator in v &&
                                    "Literal" !== t.left.type
                                ) {
                                    var l = s(t.left),
                                        c = s(t.right);
                                    a(
                                        t,
                                        "__$__(" +
                                            l +
                                            ', "' +
                                            t.operator +
                                            '", ' +
                                            c +
                                            ")"
                                    );
                                }
                                break;
                            case "UpdateExpression":
                            case "AssignmentExpression":
                                var d = e && e.type;
                                if (
                                    !(
                                        "ForStatement" === d ||
                                        ("BinaryExpression" === d &&
                                            /^[=!<>]/.test(e.operator)) ||
                                        ("MemberExpression" === d && e.computed)
                                    )
                                )
                                    if ("UpdateExpression" === t.type) {
                                        var u = s(t.argument),
                                            f =
                                                u +
                                                " = __$__(" +
                                                u +
                                                ', "' +
                                                t.operator[0] +
                                                '", 1)';
                                        t.prefix ||
                                            ("AssignmentExpression" !== d &&
                                                "VariableDeclarator" !== d) ||
                                            (f = u + "; " + f),
                                            a(t, f);
                                    } else if (
                                        /^.=$/.test(t.operator) &&
                                        "Literal" !== t.left.type
                                    ) {
                                        var l = s(t.left),
                                            c = s(t.right);
                                        a(
                                            t,
                                            l +
                                                " = __$__(" +
                                                l +
                                                ', "' +
                                                t.operator[0] +
                                                '", ' +
                                                c +
                                                ")"
                                        );
                                    }
                        }
                    }
                }
                if (!t) return "";
                (n = n || {}), (e = e || "");
                var h = [],
                    u = null,
                    l = p.version,
                    c = /\r\n|\n|\r/gm;
                if (
                    (p.chrome && l >= 30) ||
                    (p.webkit && l >= 537.76) ||
                    (p.firefox && l >= 23)
                ) {
                    var d = 0;
                    if (0 === window.location.href.indexOf(e)) {
                        var f =
                            document.getElementsByTagName("html")[0].innerHTML;
                        d = f.substr(0, f.indexOf(t) + 1).match(c).length + 1;
                    }
                    var _ = ["AAAA"];
                    (_.length = (t.match(c) || []).length + 1 + d),
                        (u = {
                            version: 3,
                            file: e,
                            names: [],
                            mappings: _.join(";AACA"),
                            sourceRoot: "",
                            sources: [e]
                        });
                    var g = n.source || (!e && t);
                    g && (u.sourcesContent = [g]);
                }
                return (
                    o(
                        i(t, {
                            ranges: !0
                        })
                    ),
                    u &&
                        (t =
                            Array(d + 1).join("\n") +
                            t +
                            "\n//# sourceMappingURL=data:application/json;base64," +
                            btoa(
                                unescape(encodeURIComponent(JSON.stringify(u)))
                            ) +
                            "\n//# sourceURL=" +
                            (e || "paperscript")),
                    t
                );
            }
            function a(i, r, a, o) {
                function u(t, e) {
                    for (var n in t)
                        (!e && /^_/.test(n)) ||
                            !RegExp(
                                "([\\b\\s\\W]|^)" +
                                    n.replace(/\$/g, "\\$") +
                                    "\\b"
                            ).test(i) ||
                            (g.push(n), v.push(t[n]));
                }
                paper = r;
                var l,
                    c = r.getView(),
                    d = /\s+on(?:Key|Mouse)(?:Up|Down|Move|Drag)\b/.test(i)
                        ? new Y()
                        : null,
                    f = d ? d._events : [],
                    _ = ["onFrame", "onResize"].concat(f),
                    g = [],
                    v = [];
                if (
                    ((i = s(i, a, o)),
                    u(
                        {
                            __$__: t,
                            $__: n,
                            paper: r,
                            view: c,
                            tool: d
                        },
                        !0
                    ),
                    u(r),
                    (_ = e
                        .each(
                            _,
                            function (t) {
                                RegExp("\\s+" + t + "\\b").test(i) &&
                                    (g.push(t), this.push(t + ": " + t));
                            },
                            []
                        )
                        .join(", ")),
                    _ && (i += "\nreturn { " + _ + " };"),
                    p.chrome || p.firefox)
                ) {
                    var m = document.createElement("script"),
                        y = document.head;
                    p.firefox && (i = "\n" + i),
                        m.appendChild(
                            document.createTextNode(
                                "paper._execute = function(" +
                                    g +
                                    ") {" +
                                    i +
                                    "\n}"
                            )
                        ),
                        y.appendChild(m),
                        (l = paper._execute),
                        delete paper._execute,
                        y.removeChild(m);
                } else l = Function(g, i);
                var w = l.apply(r, v) || {};
                e.each(f, function (t) {
                    var e = w[t];
                    e && (d[t] = e);
                }),
                    c &&
                        (w.onResize && c.setOnResize(w.onResize),
                        c.fire("resize", {
                            size: c.size,
                            delta: new h()
                        }),
                        w.onFrame && c.setOnFrame(w.onFrame),
                        c.update());
            }
            function o(t) {
                if (
                    /^text\/(?:x-|)paperscript$/.test(t.type) &&
                    "true" !== r.getAttribute(t, "ignore")
                ) {
                    var e = r.getAttribute(t, "canvas"),
                        n = document.getElementById(e),
                        i = t.src,
                        s = "data-paper-scope";
                    if (!n)
                        throw Error(
                            'Unable to find canvas with id "' + e + '"'
                        );
                    var o = r.get(n.getAttribute(s)) || new r().setup(n);
                    return (
                        n.setAttribute(s, o._id),
                        i
                            ? K.request("get", i, function (t) {
                                  a(t, o, i);
                              })
                            : a(t.innerHTML, o, t.baseURI),
                        t.setAttribute("data-paper-ignore", "true"),
                        o
                    );
                }
            }
            function u() {
                e.each(document.getElementsByTagName("script"), o);
            }
            function l(t) {
                return t ? o(t) : u();
            }
            var d,
                f,
                _ = this;
            !(function (t, e) {
                return "object" == typeof d && "object" == typeof module
                    ? e(d)
                    : "function" == typeof f && f.amd
                    ? f(["exports"], e)
                    : (e(t.acorn || (t.acorn = {})), void 0);
            })(this, function (t) {
                "use strict";
                function e(t) {
                    ce = t || {};
                    for (var e in ge)
                        Object.prototype.hasOwnProperty.call(ce, e) ||
                            (ce[e] = ge[e]);
                    _e = ce.sourceFile || null;
                }
                function n(t, e) {
                    var n = pe(de, t);
                    e += " (" + n.line + ":" + n.column + ")";
                    var i = new SyntaxError(e);
                    throw ((i.pos = t), (i.loc = n), (i.raisedAt = ve), i);
                }
                function i(t) {
                    function e(t) {
                        if (1 == t.length)
                            return (n +=
                                "return str === " + JSON.stringify(t[0]) + ";");
                        n += "switch(str){";
                        for (var e = 0; e < t.length; ++e)
                            n += "case " + JSON.stringify(t[e]) + ":";
                        n += "return true}return false;";
                    }
                    t = t.split(" ");
                    var n = "",
                        i = [];
                    t: for (var r = 0; r < t.length; ++r) {
                        for (var s = 0; s < i.length; ++s)
                            if (i[s][0].length == t[r].length) {
                                i[s].push(t[r]);
                                continue t;
                            }
                        i.push([t[r]]);
                    }
                    if (i.length > 3) {
                        i.sort(function (t, e) {
                            return e.length - t.length;
                        }),
                            (n += "switch(str.length){");
                        for (var r = 0; r < i.length; ++r) {
                            var a = i[r];
                            (n += "case " + a[0].length + ":"), e(a);
                        }
                        n += "}";
                    } else e(t);
                    return Function("str", n);
                }
                function r() {
                    (this.line = Pe), (this.column = ve - ke);
                }
                function s() {
                    (Pe = 1), (ve = ke = 0), (Se = !0), u();
                }
                function a(t, e) {
                    (ye = ve),
                        ce.locations && (xe = new r()),
                        (be = t),
                        u(),
                        (Ce = e),
                        (Se = t.beforeExpr);
                }
                function o() {
                    var t = ce.onComment && ce.locations && new r(),
                        e = ve,
                        i = de.indexOf("*/", (ve += 2));
                    if (
                        (-1 === i && n(ve - 2, "Unterminated comment"),
                        (ve = i + 2),
                        ce.locations)
                    ) {
                        Yn.lastIndex = e;
                        for (var s; (s = Yn.exec(de)) && s.index < ve; )
                            ++Pe, (ke = s.index + s[0].length);
                    }
                    ce.onComment &&
                        ce.onComment(
                            !0,
                            de.slice(e + 2, i),
                            e,
                            ve,
                            t,
                            ce.locations && new r()
                        );
                }
                function h() {
                    for (
                        var t = ve,
                            e = ce.onComment && ce.locations && new r(),
                            n = de.charCodeAt((ve += 2));
                        fe > ve &&
                        10 !== n &&
                        13 !== n &&
                        8232 !== n &&
                        8329 !== n;

                    )
                        ++ve, (n = de.charCodeAt(ve));
                    ce.onComment &&
                        ce.onComment(
                            !1,
                            de.slice(t + 2, ve),
                            t,
                            ve,
                            e,
                            ce.locations && new r()
                        );
                }
                function u() {
                    for (; fe > ve; ) {
                        var t = de.charCodeAt(ve);
                        if (32 === t) ++ve;
                        else if (13 === t) {
                            ++ve;
                            var e = de.charCodeAt(ve);
                            10 === e && ++ve, ce.locations && (++Pe, (ke = ve));
                        } else if (10 === t) ++ve, ++Pe, (ke = ve);
                        else if (14 > t && t > 8) ++ve;
                        else if (47 === t) {
                            var e = de.charCodeAt(ve + 1);
                            if (42 === e) o();
                            else {
                                if (47 !== e) break;
                                h();
                            }
                        } else if (160 === t) ++ve;
                        else {
                            if (!(t >= 5760 && Hn.test(String.fromCharCode(t))))
                                break;
                            ++ve;
                        }
                    }
                }
                function l() {
                    var t = de.charCodeAt(ve + 1);
                    return t >= 48 && 57 >= t ? S(!0) : (++ve, a(xn));
                }
                function c() {
                    var t = de.charCodeAt(ve + 1);
                    return Se ? (++ve, x()) : 61 === t ? w(Pn, 2) : w(Cn, 1);
                }
                function d() {
                    var t = de.charCodeAt(ve + 1);
                    return 61 === t ? w(Pn, 2) : w(Dn, 1);
                }
                function f(t) {
                    var e = de.charCodeAt(ve + 1);
                    return e === t
                        ? w(124 === t ? An : In, 2)
                        : 61 === e
                        ? w(Pn, 2)
                        : w(124 === t ? On : Ln, 1);
                }
                function _() {
                    var t = de.charCodeAt(ve + 1);
                    return 61 === t ? w(Pn, 2) : w(Tn, 1);
                }
                function g(t) {
                    var e = de.charCodeAt(ve + 1);
                    return e === t ? w(Mn, 2) : 61 === e ? w(Pn, 2) : w(kn, 1);
                }
                function p(t) {
                    var e = de.charCodeAt(ve + 1),
                        n = 1;
                    return e === t
                        ? ((n =
                              62 === t && 62 === de.charCodeAt(ve + 2) ? 3 : 2),
                          61 === de.charCodeAt(ve + n)
                              ? w(Pn, n + 1)
                              : w(jn, n))
                        : (61 === e &&
                              (n = 61 === de.charCodeAt(ve + 2) ? 3 : 2),
                          w(Nn, n));
                }
                function v(t) {
                    var e = de.charCodeAt(ve + 1);
                    return 61 === e
                        ? w(En, 61 === de.charCodeAt(ve + 2) ? 3 : 2)
                        : w(61 === t ? Sn : zn, 1);
                }
                function m(t) {
                    switch (t) {
                        case 46:
                            return l();
                        case 40:
                            return ++ve, a(pn);
                        case 41:
                            return ++ve, a(vn);
                        case 59:
                            return ++ve, a(yn);
                        case 44:
                            return ++ve, a(mn);
                        case 91:
                            return ++ve, a(dn);
                        case 93:
                            return ++ve, a(fn);
                        case 123:
                            return ++ve, a(_n);
                        case 125:
                            return ++ve, a(gn);
                        case 58:
                            return ++ve, a(wn);
                        case 63:
                            return ++ve, a(bn);
                        case 48:
                            var e = de.charCodeAt(ve + 1);
                            if (120 === e || 88 === e) return C();
                        case 49:
                        case 50:
                        case 51:
                        case 52:
                        case 53:
                        case 54:
                        case 55:
                        case 56:
                        case 57:
                            return S(!1);
                        case 34:
                        case 39:
                            return P(t);
                        case 47:
                            return c(t);
                        case 37:
                        case 42:
                            return d();
                        case 124:
                        case 38:
                            return f(t);
                        case 94:
                            return _();
                        case 43:
                        case 45:
                            return g(t);
                        case 60:
                        case 62:
                            return p(t);
                        case 61:
                        case 33:
                            return v(t);
                        case 126:
                            return w(zn, 1);
                    }
                    return !1;
                }
                function y(t) {
                    if (
                        (t ? (ve = me + 1) : (me = ve),
                        ce.locations && (we = new r()),
                        t)
                    )
                        return x();
                    if (ve >= fe) return a(Be);
                    var e = de.charCodeAt(ve);
                    if (Kn(e) || 92 === e) return z();
                    var i = m(e);
                    if (i === !1) {
                        var s = String.fromCharCode(e);
                        if ("\\" === s || $n.test(s)) return z();
                        n(ve, "Unexpected character '" + s + "'");
                    }
                    return i;
                }
                function w(t, e) {
                    var n = de.slice(ve, ve + e);
                    (ve += e), a(t, n);
                }
                function x() {
                    for (var t, e, i = "", r = ve; ; ) {
                        ve >= fe && n(r, "Unterminated regular expression");
                        var s = de.charAt(ve);
                        if (
                            (Jn.test(s) &&
                                n(r, "Unterminated regular expression"),
                            t)
                        )
                            t = !1;
                        else {
                            if ("[" === s) e = !0;
                            else if ("]" === s && e) e = !1;
                            else if ("/" === s && !e) break;
                            t = "\\" === s;
                        }
                        ++ve;
                    }
                    var i = de.slice(r, ve);
                    ++ve;
                    var o = M();
                    return (
                        o &&
                            !/^[gmsiy]*$/.test(o) &&
                            n(r, "Invalid regexp flag"),
                        a(Ne, RegExp(i, o))
                    );
                }
                function b(t, e) {
                    for (
                        var n = ve, i = 0, r = 0, s = null == e ? 1 / 0 : e;
                        s > r;
                        ++r
                    ) {
                        var a,
                            o = de.charCodeAt(ve);
                        if (
                            ((a =
                                o >= 97
                                    ? o - 97 + 10
                                    : o >= 65
                                    ? o - 65 + 10
                                    : o >= 48 && 57 >= o
                                    ? o - 48
                                    : 1 / 0),
                            a >= t)
                        )
                            break;
                        ++ve, (i = i * t + a);
                    }
                    return ve === n || (null != e && ve - n !== e) ? null : i;
                }
                function C() {
                    ve += 2;
                    var t = b(16);
                    return (
                        null == t && n(me + 2, "Expected hexadecimal number"),
                        Kn(de.charCodeAt(ve)) &&
                            n(ve, "Identifier directly after number"),
                        a(Ee, t)
                    );
                }
                function S(t) {
                    var e = ve,
                        i = !1,
                        r = 48 === de.charCodeAt(ve);
                    t || null !== b(10) || n(e, "Invalid number"),
                        46 === de.charCodeAt(ve) && (++ve, b(10), (i = !0));
                    var s = de.charCodeAt(ve);
                    (69 === s || 101 === s) &&
                        ((s = de.charCodeAt(++ve)),
                        (43 === s || 45 === s) && ++ve,
                        null === b(10) && n(e, "Invalid number"),
                        (i = !0)),
                        Kn(de.charCodeAt(ve)) &&
                            n(ve, "Identifier directly after number");
                    var o,
                        h = de.slice(e, ve);
                    return (
                        i
                            ? (o = parseFloat(h))
                            : r && 1 !== h.length
                            ? /[89]/.test(h) || Te
                                ? n(e, "Invalid number")
                                : (o = parseInt(h, 8))
                            : (o = parseInt(h, 10)),
                        a(Ee, o)
                    );
                }
                function P(t) {
                    ve++;
                    for (var e = ""; ; ) {
                        ve >= fe && n(me, "Unterminated string constant");
                        var i = de.charCodeAt(ve);
                        if (i === t) return ++ve, a(je, e);
                        if (92 === i) {
                            i = de.charCodeAt(++ve);
                            var r = /^[0-7]+/.exec(de.slice(ve, ve + 3));
                            for (r && (r = r[0]); r && parseInt(r, 8) > 255; )
                                r = r.slice(0, r.length - 1);
                            if (("0" === r && (r = null), ++ve, r))
                                Te && n(ve - 2, "Octal literal in strict mode"),
                                    (e += String.fromCharCode(parseInt(r, 8))),
                                    (ve += r.length - 1);
                            else
                                switch (i) {
                                    case 110:
                                        e += "\n";
                                        break;
                                    case 114:
                                        e += "\r";
                                        break;
                                    case 120:
                                        e += String.fromCharCode(k(2));
                                        break;
                                    case 117:
                                        e += String.fromCharCode(k(4));
                                        break;
                                    case 85:
                                        e += String.fromCharCode(k(8));
                                        break;
                                    case 116:
                                        e += "	";
                                        break;
                                    case 98:
                                        e += "\b";
                                        break;
                                    case 118:
                                        e += "";
                                        break;
                                    case 102:
                                        e += "\f";
                                        break;
                                    case 48:
                                        e += "\0";
                                        break;
                                    case 13:
                                        10 === de.charCodeAt(ve) && ++ve;
                                    case 10:
                                        ce.locations && ((ke = ve), ++Pe);
                                        break;
                                    default:
                                        e += String.fromCharCode(i);
                                }
                        } else
                            (13 === i ||
                                10 === i ||
                                8232 === i ||
                                8329 === i) &&
                                n(me, "Unterminated string constant"),
                                (e += String.fromCharCode(i)),
                                ++ve;
                    }
                }
                function k(t) {
                    var e = b(16, t);
                    return (
                        null === e && n(me, "Bad character escape sequence"), e
                    );
                }
                function M() {
                    Rn = !1;
                    for (var t, e = !0, i = ve; ; ) {
                        var r = de.charCodeAt(ve);
                        if (Qn(r)) Rn && (t += de.charAt(ve)), ++ve;
                        else {
                            if (92 !== r) break;
                            Rn || (t = de.slice(i, ve)),
                                (Rn = !0),
                                117 != de.charCodeAt(++ve) &&
                                    n(
                                        ve,
                                        "Expecting Unicode escape sequence \\uXXXX"
                                    ),
                                ++ve;
                            var s = k(4),
                                a = String.fromCharCode(s);
                            a || n(ve - 1, "Invalid Unicode escape"),
                                (e ? Kn(s) : Qn(s)) ||
                                    n(ve - 4, "Invalid Unicode escape"),
                                (t += a);
                        }
                        e = !1;
                    }
                    return Rn ? t : de.slice(i, ve);
                }
                function z() {
                    var t = M(),
                        e = De;
                    return (
                        Rn ||
                            (Un(t)
                                ? (e = cn[t])
                                : ((ce.forbidReserved &&
                                      (3 === ce.ecmaVersion ? Fn : Vn)(t)) ||
                                      (Te && qn(t))) &&
                                  n(me, "The keyword '" + t + "' is reserved")),
                        a(e, t)
                    );
                }
                function A() {
                    (Me = me), (ze = ye), (Ae = xe), y();
                }
                function I(t) {
                    for (Te = t, ve = ze; ke > ve; )
                        (ke = de.lastIndexOf("\n", ke - 2) + 1), --Pe;
                    u(), y();
                }
                function O() {
                    (this.type = null), (this.start = me), (this.end = null);
                }
                function T() {
                    (this.start = we),
                        (this.end = null),
                        null !== _e && (this.source = _e);
                }
                function L() {
                    var t = new O();
                    return (
                        ce.locations && (t.loc = new T()),
                        ce.ranges && (t.range = [me, 0]),
                        t
                    );
                }
                function E(t) {
                    var e = new O();
                    return (
                        (e.start = t.start),
                        ce.locations &&
                            ((e.loc = new T()), (e.loc.start = t.loc.start)),
                        ce.ranges && (e.range = [t.range[0], 0]),
                        e
                    );
                }
                function N(t, e) {
                    return (
                        (t.type = e),
                        (t.end = ze),
                        ce.locations && (t.loc.end = Ae),
                        ce.ranges && (t.range[1] = ze),
                        t
                    );
                }
                function j(t) {
                    return (
                        ce.ecmaVersion >= 5 &&
                        "ExpressionStatement" === t.type &&
                        "Literal" === t.expression.type &&
                        "use strict" === t.expression.value
                    );
                }
                function D(t) {
                    return be === t ? (A(), !0) : void 0;
                }
                function B() {
                    return (
                        !ce.strictSemicolons &&
                        (be === Be || be === gn || Jn.test(de.slice(ze, me)))
                    );
                }
                function R() {
                    D(yn) || B() || V();
                }
                function F(t) {
                    be === t ? A() : V();
                }
                function V() {
                    n(me, "Unexpected token");
                }
                function q(t) {
                    "Identifier" !== t.type &&
                        "MemberExpression" !== t.type &&
                        n(t.start, "Assigning to rvalue"),
                        Te &&
                            "Identifier" === t.type &&
                            Zn(t.name) &&
                            n(
                                t.start,
                                "Assigning to " + t.name + " in strict mode"
                            );
                }
                function Z(t) {
                    (Me = ze = ve),
                        ce.locations && (Ae = new r()),
                        (Ie = Te = null),
                        (Oe = []),
                        y();
                    var e = t || L(),
                        n = !0;
                    for (t || (e.body = []); be !== Be; ) {
                        var i = U();
                        e.body.push(i), n && j(i) && I(!0), (n = !1);
                    }
                    return N(e, "Program");
                }
                function U() {
                    be === Cn && y(!0);
                    var t = be,
                        e = L();
                    switch (t) {
                        case Re:
                        case qe:
                            A();
                            var i = t === Re;
                            D(yn) || B()
                                ? (e.label = null)
                                : be !== De
                                ? V()
                                : ((e.label = le()), R());
                            for (var r = 0; r < Oe.length; ++r) {
                                var s = Oe[r];
                                if (
                                    null == e.label ||
                                    s.name === e.label.name
                                ) {
                                    if (
                                        null != s.kind &&
                                        (i || "loop" === s.kind)
                                    )
                                        break;
                                    if (e.label && i) break;
                                }
                            }
                            return (
                                r === Oe.length &&
                                    n(e.start, "Unsyntactic " + t.keyword),
                                N(e, i ? "BreakStatement" : "ContinueStatement")
                            );
                        case Ze:
                            return A(), R(), N(e, "DebuggerStatement");
                        case He:
                            return (
                                A(),
                                Oe.push(ti),
                                (e.body = U()),
                                Oe.pop(),
                                F(nn),
                                (e.test = H()),
                                R(),
                                N(e, "DoWhileStatement")
                            );
                        case $e:
                            if ((A(), Oe.push(ti), F(pn), be === yn))
                                return G(e, null);
                            if (be === en) {
                                var a = L();
                                return (
                                    A(),
                                    X(a, !0),
                                    1 === a.declarations.length && D(ln)
                                        ? $(e, a)
                                        : G(e, a)
                                );
                            }
                            var a = J(!1, !0);
                            return D(ln) ? (q(a), $(e, a)) : G(e, a);
                        case Xe:
                            return A(), he(e, !0);
                        case Je:
                            return (
                                A(),
                                (e.test = H()),
                                (e.consequent = U()),
                                (e.alternate = D(We) ? U() : null),
                                N(e, "IfStatement")
                            );
                        case Ye:
                            return (
                                Ie || n(me, "'return' outside of function"),
                                A(),
                                D(yn) || B()
                                    ? (e.argument = null)
                                    : ((e.argument = J()), R()),
                                N(e, "ReturnStatement")
                            );
                        case Ke:
                            A(),
                                (e.discriminant = H()),
                                (e.cases = []),
                                F(_n),
                                Oe.push(ei);
                            for (var o, h; be != gn; )
                                if (be === Fe || be === Ue) {
                                    var u = be === Fe;
                                    o && N(o, "SwitchCase"),
                                        e.cases.push((o = L())),
                                        (o.consequent = []),
                                        A(),
                                        u
                                            ? (o.test = J())
                                            : (h &&
                                                  n(
                                                      Me,
                                                      "Multiple default clauses"
                                                  ),
                                              (h = !0),
                                              (o.test = null)),
                                        F(wn);
                                } else o || V(), o.consequent.push(U());
                            return (
                                o && N(o, "SwitchCase"),
                                A(),
                                Oe.pop(),
                                N(e, "SwitchStatement")
                            );
                        case Qe:
                            return (
                                A(),
                                Jn.test(de.slice(ze, me)) &&
                                    n(ze, "Illegal newline after throw"),
                                (e.argument = J()),
                                R(),
                                N(e, "ThrowStatement")
                            );
                        case tn:
                            if (
                                (A(),
                                (e.block = W()),
                                (e.handler = null),
                                be === Ve)
                            ) {
                                var l = L();
                                A(),
                                    F(pn),
                                    (l.param = le()),
                                    Te &&
                                        Zn(l.param.name) &&
                                        n(
                                            l.param.start,
                                            "Binding " +
                                                l.param.name +
                                                " in strict mode"
                                        ),
                                    F(vn),
                                    (l.guard = null),
                                    (l.body = W()),
                                    (e.handler = N(l, "CatchClause"));
                            }
                            return (
                                (e.guardedHandlers = Le),
                                (e.finalizer = D(Ge) ? W() : null),
                                e.handler ||
                                    e.finalizer ||
                                    n(
                                        e.start,
                                        "Missing catch or finally clause"
                                    ),
                                N(e, "TryStatement")
                            );
                        case en:
                            return A(), (e = X(e)), R(), e;
                        case nn:
                            return (
                                A(),
                                (e.test = H()),
                                Oe.push(ti),
                                (e.body = U()),
                                Oe.pop(),
                                N(e, "WhileStatement")
                            );
                        case rn:
                            return (
                                Te && n(me, "'with' in strict mode"),
                                A(),
                                (e.object = H()),
                                (e.body = U()),
                                N(e, "WithStatement")
                            );
                        case _n:
                            return W();
                        case yn:
                            return A(), N(e, "EmptyStatement");
                        default:
                            var c = Ce,
                                d = J();
                            if (t === De && "Identifier" === d.type && D(wn)) {
                                for (var r = 0; r < Oe.length; ++r)
                                    Oe[r].name === c &&
                                        n(
                                            d.start,
                                            "Label '" +
                                                c +
                                                "' is already declared"
                                        );
                                var f = be.isLoop
                                    ? "loop"
                                    : be === Ke
                                    ? "switch"
                                    : null;
                                return (
                                    Oe.push({
                                        name: c,
                                        kind: f
                                    }),
                                    (e.body = U()),
                                    Oe.pop(),
                                    (e.label = d),
                                    N(e, "LabeledStatement")
                                );
                            }
                            return (
                                (e.expression = d),
                                R(),
                                N(e, "ExpressionStatement")
                            );
                    }
                }
                function H() {
                    F(pn);
                    var t = J();
                    return F(vn), t;
                }
                function W(t) {
                    var e,
                        n = L(),
                        i = !0,
                        r = !1;
                    for (n.body = [], F(_n); !D(gn); ) {
                        var s = U();
                        n.body.push(s),
                            i && t && j(s) && ((e = r), I((r = !0))),
                            (i = !1);
                    }
                    return r && !e && I(!1), N(n, "BlockStatement");
                }
                function G(t, e) {
                    return (
                        (t.init = e),
                        F(yn),
                        (t.test = be === yn ? null : J()),
                        F(yn),
                        (t.update = be === vn ? null : J()),
                        F(vn),
                        (t.body = U()),
                        Oe.pop(),
                        N(t, "ForStatement")
                    );
                }
                function $(t, e) {
                    return (
                        (t.left = e),
                        (t.right = J()),
                        F(vn),
                        (t.body = U()),
                        Oe.pop(),
                        N(t, "ForInStatement")
                    );
                }
                function X(t, e) {
                    for (t.declarations = [], t.kind = "var"; ; ) {
                        var i = L();
                        if (
                            ((i.id = le()),
                            Te &&
                                Zn(i.id.name) &&
                                n(
                                    i.id.start,
                                    "Binding " + i.id.name + " in strict mode"
                                ),
                            (i.init = D(Sn) ? J(!0, e) : null),
                            t.declarations.push(N(i, "VariableDeclarator")),
                            !D(mn))
                        )
                            break;
                    }
                    return N(t, "VariableDeclaration");
                }
                function J(t, e) {
                    var n = Y(e);
                    if (!t && be === mn) {
                        var i = E(n);
                        for (i.expressions = [n]; D(mn); )
                            i.expressions.push(Y(e));
                        return N(i, "SequenceExpression");
                    }
                    return n;
                }
                function Y(t) {
                    var e = K(t);
                    if (be.isAssign) {
                        var n = E(e);
                        return (
                            (n.operator = Ce),
                            (n.left = e),
                            A(),
                            (n.right = Y(t)),
                            q(e),
                            N(n, "AssignmentExpression")
                        );
                    }
                    return e;
                }
                function K(t) {
                    var e = Q(t);
                    if (D(bn)) {
                        var n = E(e);
                        return (
                            (n.test = e),
                            (n.consequent = J(!0)),
                            F(wn),
                            (n.alternate = J(!0, t)),
                            N(n, "ConditionalExpression")
                        );
                    }
                    return e;
                }
                function Q(t) {
                    return te(ee(), -1, t);
                }
                function te(t, e, n) {
                    var i = be.binop;
                    if (null != i && (!n || be !== ln) && i > e) {
                        var r = E(t);
                        (r.left = t),
                            (r.operator = Ce),
                            A(),
                            (r.right = te(ee(), i, n));
                        var r = N(
                            r,
                            /&&|\|\|/.test(r.operator)
                                ? "LogicalExpression"
                                : "BinaryExpression"
                        );
                        return te(r, e, n);
                    }
                    return t;
                }
                function ee() {
                    if (be.prefix) {
                        var t = L(),
                            e = be.isUpdate;
                        return (
                            (t.operator = Ce),
                            (t.prefix = !0),
                            A(),
                            (t.argument = ee()),
                            e
                                ? q(t.argument)
                                : Te &&
                                  "delete" === t.operator &&
                                  "Identifier" === t.argument.type &&
                                  n(
                                      t.start,
                                      "Deleting local variable in strict mode"
                                  ),
                            N(t, e ? "UpdateExpression" : "UnaryExpression")
                        );
                    }
                    for (var i = ne(); be.postfix && !B(); ) {
                        var t = E(i);
                        (t.operator = Ce),
                            (t.prefix = !1),
                            (t.argument = i),
                            q(i),
                            A(),
                            (i = N(t, "UpdateExpression"));
                    }
                    return i;
                }
                function ne() {
                    return ie(re());
                }
                function ie(t, e) {
                    if (D(xn)) {
                        var n = E(t);
                        return (
                            (n.object = t),
                            (n.property = le(!0)),
                            (n.computed = !1),
                            ie(N(n, "MemberExpression"), e)
                        );
                    }
                    if (D(dn)) {
                        var n = E(t);
                        return (
                            (n.object = t),
                            (n.property = J()),
                            (n.computed = !0),
                            F(fn),
                            ie(N(n, "MemberExpression"), e)
                        );
                    }
                    if (!e && D(pn)) {
                        var n = E(t);
                        return (
                            (n.callee = t),
                            (n.arguments = ue(vn, !1)),
                            ie(N(n, "CallExpression"), e)
                        );
                    }
                    return t;
                }
                function re() {
                    switch (be) {
                        case an:
                            var t = L();
                            return A(), N(t, "ThisExpression");
                        case De:
                            return le();
                        case Ee:
                        case je:
                        case Ne:
                            var t = L();
                            return (
                                (t.value = Ce),
                                (t.raw = de.slice(me, ye)),
                                A(),
                                N(t, "Literal")
                            );
                        case on:
                        case hn:
                        case un:
                            var t = L();
                            return (
                                (t.value = be.atomValue),
                                (t.raw = be.keyword),
                                A(),
                                N(t, "Literal")
                            );
                        case pn:
                            var e = we,
                                n = me;
                            A();
                            var i = J();
                            return (
                                (i.start = n),
                                (i.end = ye),
                                ce.locations &&
                                    ((i.loc.start = e), (i.loc.end = xe)),
                                ce.ranges && (i.range = [n, ye]),
                                F(vn),
                                i
                            );
                        case dn:
                            var t = L();
                            return (
                                A(),
                                (t.elements = ue(fn, !0, !0)),
                                N(t, "ArrayExpression")
                            );
                        case _n:
                            return ae();
                        case Xe:
                            var t = L();
                            return A(), he(t, !1);
                        case sn:
                            return se();
                        default:
                            V();
                    }
                }
                function se() {
                    var t = L();
                    return (
                        A(),
                        (t.callee = ie(re(), !0)),
                        (t.arguments = D(pn) ? ue(vn, !1) : Le),
                        N(t, "NewExpression")
                    );
                }
                function ae() {
                    var t = L(),
                        e = !0,
                        i = !1;
                    for (t.properties = [], A(); !D(gn); ) {
                        if (e) e = !1;
                        else if ((F(mn), ce.allowTrailingCommas && D(gn)))
                            break;
                        var r,
                            s = {
                                key: oe()
                            },
                            a = !1;
                        if (
                            (D(wn)
                                ? ((s.value = J(!0)), (r = s.kind = "init"))
                                : ce.ecmaVersion >= 5 &&
                                  "Identifier" === s.key.type &&
                                  ("get" === s.key.name || "set" === s.key.name)
                                ? ((a = i = !0),
                                  (r = s.kind = s.key.name),
                                  (s.key = oe()),
                                  be !== pn && V(),
                                  (s.value = he(L(), !1)))
                                : V(),
                            "Identifier" === s.key.type && (Te || i))
                        )
                            for (var o = 0; o < t.properties.length; ++o) {
                                var h = t.properties[o];
                                if (h.key.name === s.key.name) {
                                    var u =
                                        r == h.kind ||
                                        (a && "init" === h.kind) ||
                                        ("init" === r &&
                                            ("get" === h.kind ||
                                                "set" === h.kind));
                                    u &&
                                        !Te &&
                                        "init" === r &&
                                        "init" === h.kind &&
                                        (u = !1),
                                        u &&
                                            n(
                                                s.key.start,
                                                "Redefinition of property"
                                            );
                                }
                            }
                        t.properties.push(s);
                    }
                    return N(t, "ObjectExpression");
                }
                function oe() {
                    return be === Ee || be === je ? re() : le(!0);
                }
                function he(t, e) {
                    be === De ? (t.id = le()) : e ? V() : (t.id = null),
                        (t.params = []);
                    var i = !0;
                    for (F(pn); !D(vn); )
                        i ? (i = !1) : F(mn), t.params.push(le());
                    var r = Ie,
                        s = Oe;
                    if (
                        ((Ie = !0),
                        (Oe = []),
                        (t.body = W(!0)),
                        (Ie = r),
                        (Oe = s),
                        Te || (t.body.body.length && j(t.body.body[0])))
                    )
                        for (var a = t.id ? -1 : 0; a < t.params.length; ++a) {
                            var o = 0 > a ? t.id : t.params[a];
                            if (
                                ((qn(o.name) || Zn(o.name)) &&
                                    n(
                                        o.start,
                                        "Defining '" +
                                            o.name +
                                            "' in strict mode"
                                    ),
                                a >= 0)
                            )
                                for (var h = 0; a > h; ++h)
                                    o.name === t.params[h].name &&
                                        n(
                                            o.start,
                                            "Argument name clash in strict mode"
                                        );
                        }
                    return N(
                        t,
                        e ? "FunctionDeclaration" : "FunctionExpression"
                    );
                }
                function ue(t, e, n) {
                    for (var i = [], r = !0; !D(t); ) {
                        if (r) r = !1;
                        else if ((F(mn), e && ce.allowTrailingCommas && D(t)))
                            break;
                        n && be === mn ? i.push(null) : i.push(J(!0));
                    }
                    return i;
                }
                function le(t) {
                    var e = L();
                    return (
                        (e.name =
                            be === De
                                ? Ce
                                : (t && !ce.forbidReserved && be.keyword) ||
                                  V()),
                        A(),
                        N(e, "Identifier")
                    );
                }
                t.version = "0.3.2";
                var ce, de, fe, _e;
                t.parse = function (t, n) {
                    return (
                        (de = t + ""),
                        (fe = de.length),
                        e(n),
                        s(),
                        Z(ce.program)
                    );
                };
                var ge = (t.defaultOptions = {
                        ecmaVersion: 5,
                        strictSemicolons: !1,
                        allowTrailingCommas: !0,
                        forbidReserved: !1,
                        locations: !1,
                        onComment: null,
                        ranges: !1,
                        program: null,
                        sourceFile: null
                    }),
                    pe = (t.getLineInfo = function (t, e) {
                        for (var n = 1, i = 0; ; ) {
                            Yn.lastIndex = i;
                            var r = Yn.exec(t);
                            if (!(r && r.index < e)) break;
                            ++n, (i = r.index + r[0].length);
                        }
                        return {
                            line: n,
                            column: e - i
                        };
                    });
                t.tokenize = function (t, n) {
                    function i(t) {
                        return (
                            y(t),
                            (r.start = me),
                            (r.end = ye),
                            (r.startLoc = we),
                            (r.endLoc = xe),
                            (r.type = be),
                            (r.value = Ce),
                            r
                        );
                    }
                    (de = t + ""), (fe = de.length), e(n), s();
                    var r = {};
                    return (
                        (i.jumpTo = function (t, e) {
                            if (((ve = t), ce.locations)) {
                                (Pe = 1), (ke = Yn.lastIndex = 0);
                                for (var n; (n = Yn.exec(de)) && n.index < t; )
                                    ++Pe, (ke = n.index + n[0].length);
                            }
                            (Se = e), u();
                        }),
                        i
                    );
                };
                var ve,
                    me,
                    ye,
                    we,
                    xe,
                    be,
                    Ce,
                    Se,
                    Pe,
                    ke,
                    Me,
                    ze,
                    Ae,
                    Ie,
                    Oe,
                    Te,
                    Le = [],
                    Ee = {
                        type: "num"
                    },
                    Ne = {
                        type: "regexp"
                    },
                    je = {
                        type: "string"
                    },
                    De = {
                        type: "name"
                    },
                    Be = {
                        type: "eof"
                    },
                    Re = {
                        keyword: "break"
                    },
                    Fe = {
                        keyword: "case",
                        beforeExpr: !0
                    },
                    Ve = {
                        keyword: "catch"
                    },
                    qe = {
                        keyword: "continue"
                    },
                    Ze = {
                        keyword: "debugger"
                    },
                    Ue = {
                        keyword: "default"
                    },
                    He = {
                        keyword: "do",
                        isLoop: !0
                    },
                    We = {
                        keyword: "else",
                        beforeExpr: !0
                    },
                    Ge = {
                        keyword: "finally"
                    },
                    $e = {
                        keyword: "for",
                        isLoop: !0
                    },
                    Xe = {
                        keyword: "function"
                    },
                    Je = {
                        keyword: "if"
                    },
                    Ye = {
                        keyword: "return",
                        beforeExpr: !0
                    },
                    Ke = {
                        keyword: "switch"
                    },
                    Qe = {
                        keyword: "throw",
                        beforeExpr: !0
                    },
                    tn = {
                        keyword: "try"
                    },
                    en = {
                        keyword: "var"
                    },
                    nn = {
                        keyword: "while",
                        isLoop: !0
                    },
                    rn = {
                        keyword: "with"
                    },
                    sn = {
                        keyword: "new",
                        beforeExpr: !0
                    },
                    an = {
                        keyword: "this"
                    },
                    on = {
                        keyword: "null",
                        atomValue: null
                    },
                    hn = {
                        keyword: "true",
                        atomValue: !0
                    },
                    un = {
                        keyword: "false",
                        atomValue: !1
                    },
                    ln = {
                        keyword: "in",
                        binop: 7,
                        beforeExpr: !0
                    },
                    cn = {
                        break: Re,
                        case: Fe,
                        catch: Ve,
                        continue: qe,
                        debugger: Ze,
                        default: Ue,
                        do: He,
                        else: We,
                        finally: Ge,
                        for: $e,
                        function: Xe,
                        if: Je,
                        return: Ye,
                        switch: Ke,
                        throw: Qe,
                        try: tn,
                        var: en,
                        while: nn,
                        with: rn,
                        null: on,
                        true: hn,
                        false: un,
                        new: sn,
                        in: ln,
                        instanceof: {
                            keyword: "instanceof",
                            binop: 7,
                            beforeExpr: !0
                        },
                        this: an,
                        typeof: {
                            keyword: "typeof",
                            prefix: !0,
                            beforeExpr: !0
                        },
                        void: {
                            keyword: "void",
                            prefix: !0,
                            beforeExpr: !0
                        },
                        delete: {
                            keyword: "delete",
                            prefix: !0,
                            beforeExpr: !0
                        }
                    },
                    dn = {
                        type: "[",
                        beforeExpr: !0
                    },
                    fn = {
                        type: "]"
                    },
                    _n = {
                        type: "{",
                        beforeExpr: !0
                    },
                    gn = {
                        type: "}"
                    },
                    pn = {
                        type: "(",
                        beforeExpr: !0
                    },
                    vn = {
                        type: ")"
                    },
                    mn = {
                        type: ",",
                        beforeExpr: !0
                    },
                    yn = {
                        type: ";",
                        beforeExpr: !0
                    },
                    wn = {
                        type: ":",
                        beforeExpr: !0
                    },
                    xn = {
                        type: "."
                    },
                    bn = {
                        type: "?",
                        beforeExpr: !0
                    },
                    Cn = {
                        binop: 10,
                        beforeExpr: !0
                    },
                    Sn = {
                        isAssign: !0,
                        beforeExpr: !0
                    },
                    Pn = {
                        isAssign: !0,
                        beforeExpr: !0
                    },
                    kn = {
                        binop: 9,
                        prefix: !0,
                        beforeExpr: !0
                    },
                    Mn = {
                        postfix: !0,
                        prefix: !0,
                        isUpdate: !0
                    },
                    zn = {
                        prefix: !0,
                        beforeExpr: !0
                    },
                    An = {
                        binop: 1,
                        beforeExpr: !0
                    },
                    In = {
                        binop: 2,
                        beforeExpr: !0
                    },
                    On = {
                        binop: 3,
                        beforeExpr: !0
                    },
                    Tn = {
                        binop: 4,
                        beforeExpr: !0
                    },
                    Ln = {
                        binop: 5,
                        beforeExpr: !0
                    },
                    En = {
                        binop: 6,
                        beforeExpr: !0
                    },
                    Nn = {
                        binop: 7,
                        beforeExpr: !0
                    },
                    jn = {
                        binop: 8,
                        beforeExpr: !0
                    },
                    Dn = {
                        binop: 10,
                        beforeExpr: !0
                    };
                t.tokTypes = {
                    bracketL: dn,
                    bracketR: fn,
                    braceL: _n,
                    braceR: gn,
                    parenL: pn,
                    parenR: vn,
                    comma: mn,
                    semi: yn,
                    colon: wn,
                    dot: xn,
                    question: bn,
                    slash: Cn,
                    eq: Sn,
                    name: De,
                    eof: Be,
                    num: Ee,
                    regexp: Ne,
                    string: je
                };
                for (var Bn in cn) t.tokTypes["_" + Bn] = cn[Bn];
                var Rn,
                    Fn = i(
                        "abstract boolean byte char class double enum export extends final float goto implements import int interface long native package private protected public short static super synchronized throws transient volatile"
                    ),
                    Vn = i("class enum extends super const export import"),
                    qn = i(
                        "implements interface let package private protected public static yield"
                    ),
                    Zn = i("eval arguments"),
                    Un = i(
                        "break case catch continue debugger default do else finally for function if return switch throw try var while with null true false instanceof typeof void delete new in this"
                    ),
                    Hn =
                        /[\u1680\u180e\u2000-\u200a\u2028\u2029\u202f\u205f\u3000\ufeff]/,
                    Wn =
                        "\xaa\xb5\xba\xc0-\xd6\xd8-\xf6\xf8-\u02c1\u02c6-\u02d1\u02e0-\u02e4\u02ec\u02ee\u0370-\u0374\u0376\u0377\u037a-\u037d\u0386\u0388-\u038a\u038c\u038e-\u03a1\u03a3-\u03f5\u03f7-\u0481\u048a-\u0527\u0531-\u0556\u0559\u0561-\u0587\u05d0-\u05ea\u05f0-\u05f2\u0620-\u064a\u066e\u066f\u0671-\u06d3\u06d5\u06e5\u06e6\u06ee\u06ef\u06fa-\u06fc\u06ff\u0710\u0712-\u072f\u074d-\u07a5\u07b1\u07ca-\u07ea\u07f4\u07f5\u07fa\u0800-\u0815\u081a\u0824\u0828\u0840-\u0858\u08a0\u08a2-\u08ac\u0904-\u0939\u093d\u0950\u0958-\u0961\u0971-\u0977\u0979-\u097f\u0985-\u098c\u098f\u0990\u0993-\u09a8\u09aa-\u09b0\u09b2\u09b6-\u09b9\u09bd\u09ce\u09dc\u09dd\u09df-\u09e1\u09f0\u09f1\u0a05-\u0a0a\u0a0f\u0a10\u0a13-\u0a28\u0a2a-\u0a30\u0a32\u0a33\u0a35\u0a36\u0a38\u0a39\u0a59-\u0a5c\u0a5e\u0a72-\u0a74\u0a85-\u0a8d\u0a8f-\u0a91\u0a93-\u0aa8\u0aaa-\u0ab0\u0ab2\u0ab3\u0ab5-\u0ab9\u0abd\u0ad0\u0ae0\u0ae1\u0b05-\u0b0c\u0b0f\u0b10\u0b13-\u0b28\u0b2a-\u0b30\u0b32\u0b33\u0b35-\u0b39\u0b3d\u0b5c\u0b5d\u0b5f-\u0b61\u0b71\u0b83\u0b85-\u0b8a\u0b8e-\u0b90\u0b92-\u0b95\u0b99\u0b9a\u0b9c\u0b9e\u0b9f\u0ba3\u0ba4\u0ba8-\u0baa\u0bae-\u0bb9\u0bd0\u0c05-\u0c0c\u0c0e-\u0c10\u0c12-\u0c28\u0c2a-\u0c33\u0c35-\u0c39\u0c3d\u0c58\u0c59\u0c60\u0c61\u0c85-\u0c8c\u0c8e-\u0c90\u0c92-\u0ca8\u0caa-\u0cb3\u0cb5-\u0cb9\u0cbd\u0cde\u0ce0\u0ce1\u0cf1\u0cf2\u0d05-\u0d0c\u0d0e-\u0d10\u0d12-\u0d3a\u0d3d\u0d4e\u0d60\u0d61\u0d7a-\u0d7f\u0d85-\u0d96\u0d9a-\u0db1\u0db3-\u0dbb\u0dbd\u0dc0-\u0dc6\u0e01-\u0e30\u0e32\u0e33\u0e40-\u0e46\u0e81\u0e82\u0e84\u0e87\u0e88\u0e8a\u0e8d\u0e94-\u0e97\u0e99-\u0e9f\u0ea1-\u0ea3\u0ea5\u0ea7\u0eaa\u0eab\u0ead-\u0eb0\u0eb2\u0eb3\u0ebd\u0ec0-\u0ec4\u0ec6\u0edc-\u0edf\u0f00\u0f40-\u0f47\u0f49-\u0f6c\u0f88-\u0f8c\u1000-\u102a\u103f\u1050-\u1055\u105a-\u105d\u1061\u1065\u1066\u106e-\u1070\u1075-\u1081\u108e\u10a0-\u10c5\u10c7\u10cd\u10d0-\u10fa\u10fc-\u1248\u124a-\u124d\u1250-\u1256\u1258\u125a-\u125d\u1260-\u1288\u128a-\u128d\u1290-\u12b0\u12b2-\u12b5\u12b8-\u12be\u12c0\u12c2-\u12c5\u12c8-\u12d6\u12d8-\u1310\u1312-\u1315\u1318-\u135a\u1380-\u138f\u13a0-\u13f4\u1401-\u166c\u166f-\u167f\u1681-\u169a\u16a0-\u16ea\u16ee-\u16f0\u1700-\u170c\u170e-\u1711\u1720-\u1731\u1740-\u1751\u1760-\u176c\u176e-\u1770\u1780-\u17b3\u17d7\u17dc\u1820-\u1877\u1880-\u18a8\u18aa\u18b0-\u18f5\u1900-\u191c\u1950-\u196d\u1970-\u1974\u1980-\u19ab\u19c1-\u19c7\u1a00-\u1a16\u1a20-\u1a54\u1aa7\u1b05-\u1b33\u1b45-\u1b4b\u1b83-\u1ba0\u1bae\u1baf\u1bba-\u1be5\u1c00-\u1c23\u1c4d-\u1c4f\u1c5a-\u1c7d\u1ce9-\u1cec\u1cee-\u1cf1\u1cf5\u1cf6\u1d00-\u1dbf\u1e00-\u1f15\u1f18-\u1f1d\u1f20-\u1f45\u1f48-\u1f4d\u1f50-\u1f57\u1f59\u1f5b\u1f5d\u1f5f-\u1f7d\u1f80-\u1fb4\u1fb6-\u1fbc\u1fbe\u1fc2-\u1fc4\u1fc6-\u1fcc\u1fd0-\u1fd3\u1fd6-\u1fdb\u1fe0-\u1fec\u1ff2-\u1ff4\u1ff6-\u1ffc\u2071\u207f\u2090-\u209c\u2102\u2107\u210a-\u2113\u2115\u2119-\u211d\u2124\u2126\u2128\u212a-\u212d\u212f-\u2139\u213c-\u213f\u2145-\u2149\u214e\u2160-\u2188\u2c00-\u2c2e\u2c30-\u2c5e\u2c60-\u2ce4\u2ceb-\u2cee\u2cf2\u2cf3\u2d00-\u2d25\u2d27\u2d2d\u2d30-\u2d67\u2d6f\u2d80-\u2d96\u2da0-\u2da6\u2da8-\u2dae\u2db0-\u2db6\u2db8-\u2dbe\u2dc0-\u2dc6\u2dc8-\u2dce\u2dd0-\u2dd6\u2dd8-\u2dde\u2e2f\u3005-\u3007\u3021-\u3029\u3031-\u3035\u3038-\u303c\u3041-\u3096\u309d-\u309f\u30a1-\u30fa\u30fc-\u30ff\u3105-\u312d\u3131-\u318e\u31a0-\u31ba\u31f0-\u31ff\u3400-\u4db5\u4e00-\u9fcc\ua000-\ua48c\ua4d0-\ua4fd\ua500-\ua60c\ua610-\ua61f\ua62a\ua62b\ua640-\ua66e\ua67f-\ua697\ua6a0-\ua6ef\ua717-\ua71f\ua722-\ua788\ua78b-\ua78e\ua790-\ua793\ua7a0-\ua7aa\ua7f8-\ua801\ua803-\ua805\ua807-\ua80a\ua80c-\ua822\ua840-\ua873\ua882-\ua8b3\ua8f2-\ua8f7\ua8fb\ua90a-\ua925\ua930-\ua946\ua960-\ua97c\ua984-\ua9b2\ua9cf\uaa00-\uaa28\uaa40-\uaa42\uaa44-\uaa4b\uaa60-\uaa76\uaa7a\uaa80-\uaaaf\uaab1\uaab5\uaab6\uaab9-\uaabd\uaac0\uaac2\uaadb-\uaadd\uaae0-\uaaea\uaaf2-\uaaf4\uab01-\uab06\uab09-\uab0e\uab11-\uab16\uab20-\uab26\uab28-\uab2e\uabc0-\uabe2\uac00-\ud7a3\ud7b0-\ud7c6\ud7cb-\ud7fb\uf900-\ufa6d\ufa70-\ufad9\ufb00-\ufb06\ufb13-\ufb17\ufb1d\ufb1f-\ufb28\ufb2a-\ufb36\ufb38-\ufb3c\ufb3e\ufb40\ufb41\ufb43\ufb44\ufb46-\ufbb1\ufbd3-\ufd3d\ufd50-\ufd8f\ufd92-\ufdc7\ufdf0-\ufdfb\ufe70-\ufe74\ufe76-\ufefc\uff21-\uff3a\uff41-\uff5a\uff66-\uffbe\uffc2-\uffc7\uffca-\uffcf\uffd2-\uffd7\uffda-\uffdc",
                    Gn =
                        "\u0300-\u036f\u0483-\u0487\u0591-\u05bd\u05bf\u05c1\u05c2\u05c4\u05c5\u05c7\u0610-\u061a\u0620-\u0649\u0672-\u06d3\u06e7-\u06e8\u06fb-\u06fc\u0730-\u074a\u0800-\u0814\u081b-\u0823\u0825-\u0827\u0829-\u082d\u0840-\u0857\u08e4-\u08fe\u0900-\u0903\u093a-\u093c\u093e-\u094f\u0951-\u0957\u0962-\u0963\u0966-\u096f\u0981-\u0983\u09bc\u09be-\u09c4\u09c7\u09c8\u09d7\u09df-\u09e0\u0a01-\u0a03\u0a3c\u0a3e-\u0a42\u0a47\u0a48\u0a4b-\u0a4d\u0a51\u0a66-\u0a71\u0a75\u0a81-\u0a83\u0abc\u0abe-\u0ac5\u0ac7-\u0ac9\u0acb-\u0acd\u0ae2-\u0ae3\u0ae6-\u0aef\u0b01-\u0b03\u0b3c\u0b3e-\u0b44\u0b47\u0b48\u0b4b-\u0b4d\u0b56\u0b57\u0b5f-\u0b60\u0b66-\u0b6f\u0b82\u0bbe-\u0bc2\u0bc6-\u0bc8\u0bca-\u0bcd\u0bd7\u0be6-\u0bef\u0c01-\u0c03\u0c46-\u0c48\u0c4a-\u0c4d\u0c55\u0c56\u0c62-\u0c63\u0c66-\u0c6f\u0c82\u0c83\u0cbc\u0cbe-\u0cc4\u0cc6-\u0cc8\u0cca-\u0ccd\u0cd5\u0cd6\u0ce2-\u0ce3\u0ce6-\u0cef\u0d02\u0d03\u0d46-\u0d48\u0d57\u0d62-\u0d63\u0d66-\u0d6f\u0d82\u0d83\u0dca\u0dcf-\u0dd4\u0dd6\u0dd8-\u0ddf\u0df2\u0df3\u0e34-\u0e3a\u0e40-\u0e45\u0e50-\u0e59\u0eb4-\u0eb9\u0ec8-\u0ecd\u0ed0-\u0ed9\u0f18\u0f19\u0f20-\u0f29\u0f35\u0f37\u0f39\u0f41-\u0f47\u0f71-\u0f84\u0f86-\u0f87\u0f8d-\u0f97\u0f99-\u0fbc\u0fc6\u1000-\u1029\u1040-\u1049\u1067-\u106d\u1071-\u1074\u1082-\u108d\u108f-\u109d\u135d-\u135f\u170e-\u1710\u1720-\u1730\u1740-\u1750\u1772\u1773\u1780-\u17b2\u17dd\u17e0-\u17e9\u180b-\u180d\u1810-\u1819\u1920-\u192b\u1930-\u193b\u1951-\u196d\u19b0-\u19c0\u19c8-\u19c9\u19d0-\u19d9\u1a00-\u1a15\u1a20-\u1a53\u1a60-\u1a7c\u1a7f-\u1a89\u1a90-\u1a99\u1b46-\u1b4b\u1b50-\u1b59\u1b6b-\u1b73\u1bb0-\u1bb9\u1be6-\u1bf3\u1c00-\u1c22\u1c40-\u1c49\u1c5b-\u1c7d\u1cd0-\u1cd2\u1d00-\u1dbe\u1e01-\u1f15\u200c\u200d\u203f\u2040\u2054\u20d0-\u20dc\u20e1\u20e5-\u20f0\u2d81-\u2d96\u2de0-\u2dff\u3021-\u3028\u3099\u309a\ua640-\ua66d\ua674-\ua67d\ua69f\ua6f0-\ua6f1\ua7f8-\ua800\ua806\ua80b\ua823-\ua827\ua880-\ua881\ua8b4-\ua8c4\ua8d0-\ua8d9\ua8f3-\ua8f7\ua900-\ua909\ua926-\ua92d\ua930-\ua945\ua980-\ua983\ua9b3-\ua9c0\uaa00-\uaa27\uaa40-\uaa41\uaa4c-\uaa4d\uaa50-\uaa59\uaa7b\uaae0-\uaae9\uaaf2-\uaaf3\uabc0-\uabe1\uabec\uabed\uabf0-\uabf9\ufb20-\ufb28\ufe00-\ufe0f\ufe20-\ufe26\ufe33\ufe34\ufe4d-\ufe4f\uff10-\uff19\uff3f",
                    $n = RegExp("[" + Wn + "]"),
                    Xn = RegExp("[" + Wn + Gn + "]"),
                    Jn = /[\n\r\u2028\u2029]/,
                    Yn = /\r\n|[\n\r\u2028\u2029]/g,
                    Kn = (t.isIdentifierStart = function (t) {
                        return 65 > t
                            ? 36 === t
                            : 91 > t
                            ? !0
                            : 97 > t
                            ? 95 === t
                            : 123 > t
                            ? !0
                            : t >= 170 && $n.test(String.fromCharCode(t));
                    }),
                    Qn = (t.isIdentifierChar = function (t) {
                        return 48 > t
                            ? 36 === t
                            : 58 > t
                            ? !0
                            : 65 > t
                            ? !1
                            : 91 > t
                            ? !0
                            : 97 > t
                            ? 95 === t
                            : 123 > t
                            ? !0
                            : t >= 170 && Xn.test(String.fromCharCode(t));
                    }),
                    ti = {
                        kind: "loop"
                    },
                    ei = {
                        kind: "switch"
                    };
            });
            var g = navigator.userAgent,
                p = {};
            g.toLowerCase().replace(
                /(opera|chrome|safari|webkit|firefox|msie|trident)\/?\s*([.\d]+)(?:.*version\/([.\d]+))?(?:.*rv\:([.\d]+))?/g,
                function (t, e, n, i, r) {
                    if (!p.chrome) {
                        var s = "opera" === e ? i : n;
                        "trident" === e && ((s = r), (e = "msie")),
                            (p.version = parseFloat(s)),
                            (p.name = e),
                            (p[e] = !0),
                            p.chrome && delete p.webkit;
                    }
                }
            );
            var v = {
                    "+": "__add",
                    "-": "__subtract",
                    "*": "__multiply",
                    "/": "__divide",
                    "%": "__modulo",
                    "==": "equals",
                    "!=": "equals"
                },
                m = {
                    "-": "__negate",
                    "+": null
                },
                y = e.each(
                    [
                        "add",
                        "subtract",
                        "multiply",
                        "divide",
                        "modulo",
                        "negate"
                    ],
                    function (t) {
                        this["__" + t] = "#" + t;
                    },
                    {}
                );
            return (
                h.inject(y),
                c.inject(y),
                D.inject(y),
                "complete" === document.readyState
                    ? setTimeout(u)
                    : q.add(window, {
                          load: u
                      }),
                {
                    compile: s,
                    execute: a,
                    load: l,
                    parse: i
                }
            );
        }.call(this)),
        (paper = new (r.inject(e.exports, {
            enumerable: !0,
            Base: e,
            Numerical: o,
            DomElement: V,
            DomEvent: q,
            Http: K,
            Key: G
        }))()),
        "function" == typeof define && define.amd
            ? define("paper", paper)
            : "object" == typeof module &&
              module &&
              "object" == typeof module.exports &&
              (module.exports = paper),
        paper
    );
})();
