using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq.Expressions;
using System.Reflection;
using System.Text;
using System.Xml;
using System.Xml.Linq;
using System.Xml.Serialization;

namespace Salsa.Core.Configuration
{
    internal static class Utility
    {
        /// <summary>
        /// Finds specified XmlElement and converts it to its .Net object representation
        /// </summary>
        /// <param name="any">List of XmlElements</param>
        /// <param name="type">.Net type to convert to</param>
        /// <param name="elementName">Name of element to find</param>
        /// <param name="defaultNamespace">Namespace of element to find</param>
        /// <returns>.Net object represenation of requested XML element</returns>
        /// <remarks>Creation of XmlSerialzier's dynamically can be slow. Consider caching serializers.</remarks>
        public static object GetElementValue(XmlElement[] any, Type type, string elementName, string defaultNamespace)
        {
            object ret = null;

            // Loop through extended elements and look for element
            foreach (XmlElement xmlElement in any)
            {
                // Check element's name (wo/ ns prefix) and namespace
                if ((0 == string.Compare(xmlElement.LocalName, elementName, true, CultureInfo.InvariantCulture)) &&
                    (0 == string.Compare(xmlElement.NamespaceURI, defaultNamespace, true, CultureInfo.InvariantCulture)))
                {
                    // If found, deserialize
                    using (var reader = new XmlTextReader(new StringReader(xmlElement.OuterXml)))
                    {
                        var serializer = new XmlSerializer(
                            type,
                            new XmlAttributeOverrides(),
                            new Type[] {},
                            new XmlRootAttribute(elementName),
                            defaultNamespace);

                        ret = serializer.Deserialize(reader);
                    }

                    break;
                }
            }

            return ret;
        }

        /// <summary>
        /// Converts .Net object to an XML string
        /// </summary>
        /// <param name="o">.Net object to convert</param>
        /// <param name="elementName">Name to give root XML element</param>
        /// <param name="defaultNamespace">Namespace for the XML infoset</param>
        /// <returns>Object's XML serialization string</returns>
        public static string ConvertToXml(this object o, string elementName, string defaultNamespace)
        {
            // Serialize .Net object to XML string
            var sb = new StringBuilder();
            var serializer = new XmlSerializer(
                o.GetType(),
                new XmlAttributeOverrides(),
                new Type[] {},
                new XmlRootAttribute(elementName),
                defaultNamespace);

            using (var writer = new StringWriter(sb, NumberFormatInfo.InvariantInfo))
            {
                serializer.Serialize(writer, o);
            }

            return sb.ToString();
        }

        /// <summary>
        /// Converts .Net object to an XML string
        /// </summary>
        /// <param name="o">.Net object to convert</param>
        /// <param name="defaultNamespace">Namespace to use when serializing</param>
        /// <returns>Object's XML serialization string</returns>
        public static string ConvertToXml(this object o, string defaultNamespace)
        {
            return ConvertToXml(o, o.GetType().Name, defaultNamespace);
        }

        /// </summary>
        /// <param name="o">.Net object to convert</param>
        /// <param name="elementName">Name to give root XML element</param>
        /// <param name="defaultNamespace">Namespace for the XML infoset</param>
        /// <returns>XmlElement containing Xml serialized object</returns>
        public static XmlElement ConvertToXmlElement(this object o, string elementName, string defaultNamespace)
        {
            string xml = ConvertToXml(o, elementName, defaultNamespace);

            // Convert XML string to XmlElement
            var document = new XmlDocument();
            document.LoadXml(xml);
            return document.DocumentElement;
        }

        /// <summary>
        /// Returns the properties of the given object as XElements.
        /// Properties with null values are still returned, but as empty
        /// elements. Underscores in property names are replaces with hyphens.
        /// </summary>
        public static IEnumerable<XElement> ConvertPropertiesToXmlElements(this object source)
        {
            foreach (PropertyInfo prop in source.GetType().GetProperties())
            {
                object value = prop.GetValue(source, null);
                yield return new XElement(prop.Name.Replace("_", "-"), value);
            }
        }

        /// <summary>
        /// Returns the properties of the given object as XElements.
        /// Properties with null values are returned as empty attributes.
        /// Underscores in property names are replaces with hyphens.
        /// </summary>
        public static IEnumerable<XAttribute> ConvertPropertiesToXmlAttributes(this object source)
        {
            foreach (PropertyInfo prop in source.GetType().GetProperties())
            {
                object value = prop.GetValue(source, null);
                yield return new XAttribute(prop.Name.Replace("_", "-"), value ?? "");
            }
        }

        /// <summary>
        /// Obtains a delegate to invoke a parameterless constructor
        /// </summary>
        /// <typeparam name="TResult">The base/interface type to yield as the
        /// new value; often object except for factory pattern implementations</typeparam>
        /// <param name="type">The Type to be created</param>
        /// <returns>A delegate to the constructor if found, else null</returns>
        public static Func<TResult> Ctor<TResult>(this Type type)
        {
            ConstructorInfo ci = GetConstructor(type, Type.EmptyTypes);
            return Expression.Lambda<Func<TResult>>(
                Expression.New(ci)).Compile();
        }

        /// <summary>
        /// Obtains a delegate to invoke a constructor which takes a parameter
        /// </summary>
        /// <typeparam name="TArg1">The type of the constructor parameter</typeparam>
        /// <typeparam name="TResult">The base/interface type to yield as the
        /// new value; often object except for factory pattern implementations</typeparam>
        /// <param name="type">The Type to be created</param>
        /// <returns>A delegate to the constructor if found, else null</returns>
        public static Func<TArg1, TResult> Ctor<TArg1, TResult>(this Type type)
        {
            ConstructorInfo ci = GetConstructor(type, typeof (TArg1));
            ParameterExpression
                param1 = Expression.Parameter(typeof (TArg1), "arg1");

            return Expression.Lambda<Func<TArg1, TResult>>(
                Expression.New(ci, param1), param1).Compile();
        }

        /// <summary>
        /// Obtains a delegate to invoke a constructor with multiple parameters
        /// </summary>
        /// <typeparam name="TArg1">The type of the first constructor parameter</typeparam>
        /// <typeparam name="TArg2">The type of the second constructor parameter</typeparam>
        /// <typeparam name="TResult">The base/interface type to yield as the
        /// new value; often object except for factory pattern implementations</typeparam>
        /// <param name="type">The Type to be created</param>
        /// <returns>A delegate to the constructor if found, else null</returns>
        public static Func<TArg1, TArg2, TResult> Ctor<TArg1, TArg2, TResult>(this Type type)
        {
            ConstructorInfo ci = GetConstructor(type, typeof (TArg1), typeof (TArg2));
            ParameterExpression
                param1 = Expression.Parameter(typeof (TArg1), "arg1"),
                param2 = Expression.Parameter(typeof (TArg2), "arg2");

            return Expression.Lambda<Func<TArg1, TArg2, TResult>>(
                Expression.New(ci, param1, param2), param1, param2).Compile();
        }

        /// <summary>
        /// Obtains a delegate to invoke a constructor with multiple parameters
        /// </summary>
        /// <typeparam name="TArg1">The type of the first constructor parameter</typeparam>
        /// <typeparam name="TArg2">The type of the second constructor parameter</typeparam>
        /// <typeparam name="TArg3">The type of the third constructor parameter</typeparam>
        /// <typeparam name="TResult">The base/interface type to yield as the
        /// new value; often object except for factory pattern implementations</typeparam>
        /// <param name="type">The Type to be created</param>
        /// <returns>A delegate to the constructor if found, else null</returns>
        public static Func<TArg1, TArg2, TArg3, TResult> Ctor<TArg1, TArg2, TArg3, TResult>(this Type type)
        {
            ConstructorInfo ci = GetConstructor(type, typeof (TArg1), typeof (TArg2), typeof (TArg3));
            ParameterExpression
                param1 = Expression.Parameter(typeof (TArg1), "arg1"),
                param2 = Expression.Parameter(typeof (TArg2), "arg2"),
                param3 = Expression.Parameter(typeof (TArg3), "arg3");

            return Expression.Lambda<Func<TArg1, TArg2, TArg3, TResult>>(
                Expression.New(ci, param1, param2, param3),
                param1, param2, param3).Compile();
        }

        /// <summary>
        /// Obtains a delegate to invoke a constructor with multiple parameters
        /// </summary>
        /// <typeparam name="TArg1">The type of the first constructor parameter</typeparam>
        /// <typeparam name="TArg2">The type of the second constructor parameter</typeparam>
        /// <typeparam name="TArg3">The type of the third constructor parameter</typeparam>
        /// <typeparam name="TArg4">The type of the fourth constructor parameter</typeparam>
        /// <typeparam name="TResult">The base/interface type to yield as the
        /// new value; often object except for factory pattern implementations</typeparam>
        /// <param name="type">The Type to be created</param>
        /// <returns>A delegate to the constructor if found, else null</returns>
        public static Func<TArg1, TArg2, TArg3, TArg4, TResult> Ctor<TArg1, TArg2, TArg3, TArg4, TResult>(this Type type)
        {
            ConstructorInfo ci = GetConstructor(type, typeof (TArg1), typeof (TArg2), typeof (TArg3), typeof (TArg4));
            ParameterExpression
                param1 = Expression.Parameter(typeof (TArg1), "arg1"),
                param2 = Expression.Parameter(typeof (TArg2), "arg2"),
                param3 = Expression.Parameter(typeof (TArg3), "arg3"),
                param4 = Expression.Parameter(typeof (TArg4), "arg4");

            return Expression.Lambda<Func<TArg1, TArg2, TArg3, TArg4, TResult>>(
                Expression.New(ci, param1, param2, param3, param4),
                param1, param2, param3, param4).Compile();
        }

        private static ConstructorInfo GetConstructor(Type type, params Type[] argumentTypes)
        {
            type.ThrowIfNull("type");
            argumentTypes.ThrowIfNull("argumentTypes");

            ConstructorInfo ci = type.GetConstructor(argumentTypes);
            if (ci == null)
            {
                var sb = new StringBuilder();
                sb.Append(type.Name).Append(" has no ctor(");
                for (int i = 0; i < argumentTypes.Length; i++)
                {
                    if (i > 0)
                    {
                        sb.Append(',');
                    }
                    sb.Append(argumentTypes[i].Name);
                }
                sb.Append(')');
                throw new InvalidOperationException(sb.ToString());
            }
            return ci;
        }

        /// <summary>
        /// Throws an ArgumentNullException if the given data item is null.
        /// </summary>
        /// <param name="data">The item to check for nullity.</param>
        /// <param name="name">The name to use when throwing an exception, if necessary</param>
        public static void ThrowIfNull<T>(this T data, string name) where T : class
        {
            if (data == null)
            {
                throw new ArgumentNullException(name);
            }
        }

        /// <summary>
        /// Throws an ArgumentNullException if the given data item is null.
        /// No parameter name is specified.
        /// </summary>
        /// <param name="data">The item to check for nullity.</param>
        public static void ThrowIfNull<T>(this T data) where T : class
        {
            if (data == null)
            {
                throw new ArgumentNullException();
            }
        }
    }
}