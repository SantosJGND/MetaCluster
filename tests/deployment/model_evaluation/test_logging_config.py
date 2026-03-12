"""
Tests for deployment/model_evaluation/logging_config.py
"""
import pytest
import logging
import json
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from deployment.model_evaluation import logging_config as lc


class TestJSONFormatter:
    """Tests for JSONFormatter class."""

    def test_format_includes_timestamp(self):
        """Output includes ISO timestamp."""
        from deployment.model_evaluation.logging_config import JSONFormatter

        formatter = JSONFormatter()
        record = logging.LogRecord(
            name='test', level=logging.INFO, pathname='test.py',
            lineno=1, msg='test message', args=(), exc_info=None
        )

        result = formatter.format(record)
        parsed = json.loads(result)

        assert 'timestamp' in parsed
        assert parsed['level'] == 'INFO'
        assert parsed['message'] == 'test message'

    def test_format_includes_level(self):
        """Output includes log level."""
        from deployment.model_evaluation.logging_config import JSONFormatter

        formatter = JSONFormatter()
        record = logging.LogRecord(
            name='test', level=logging.WARNING, pathname='test.py',
            lineno=1, msg='warning message', args=(), exc_info=None
        )

        result = formatter.format(record)
        parsed = json.loads(result)

        assert parsed['level'] == 'WARNING'

    def test_format_includes_message(self):
        """Output includes message."""
        from deployment.model_evaluation.logging_config import JSONFormatter

        formatter = JSONFormatter()
        record = logging.LogRecord(
            name='test', level=logging.INFO, pathname='test.py',
            lineno=1, msg='hello world', args=(), exc_info=None
        )

        result = formatter.format(record)
        parsed = json.loads(result)

        assert parsed['message'] == 'hello world'

    def test_format_with_extra_data(self):
        """Include extra data when present."""
        from deployment.model_evaluation.logging_config import JSONFormatter

        formatter = JSONFormatter(include_extra=True)

        record = logging.LogRecord(
            name='test', level=logging.INFO, pathname='test.py',
            lineno=1, msg='test', args=(), exc_info=None
        )
        record.extra_data = {'user': 'test_user', 'action': 'login'}

        result = formatter.format(record)
        parsed = json.loads(result)

        assert 'extra' in parsed
        assert parsed['extra']['user'] == 'test_user'

    def test_format_with_exception(self):
        """Include exception info when present."""
        from deployment.model_evaluation.logging_config import JSONFormatter

        formatter = JSONFormatter()

        try:
            raise ValueError("test error")
        except ValueError:
            import sys
            record = logging.LogRecord(
                name='test', level=logging.ERROR, pathname='test.py',
                lineno=1, msg='error occurred', args=(), exc_info=sys.exc_info()
            )

        result = formatter.format(record)
        parsed = json.loads(result)

        assert 'exception' in parsed

    def test_format_includes_module_and_function(self):
        """Include module and function name."""
        formatter = lc.JSONFormatter()
        record = logging.LogRecord(
            name='test', level=logging.INFO, pathname='test.py',
            lineno=1, msg='test message', args=(), exc_info=None
        )

        result = formatter.format(record)
        parsed = json.loads(result)

        assert 'module' in parsed
        assert 'function' in parsed


class TestStructuredLogger:
    """Tests for StructuredLogger class."""

    def test_debug_logs_message(self):
        """Debug level logging."""
        from deployment.model_evaluation.logging_config import StructuredLogger

        mock_logger = Mock()
        logger = StructuredLogger('test', mock_logger)

        logger.debug('debug message')

        mock_logger.log.assert_called_once()
        call_args = mock_logger.log.call_args
        assert call_args[0][0] == logging.DEBUG

    def test_info_logs_message(self):
        """Info level logging."""
        from deployment.model_evaluation.logging_config import StructuredLogger

        mock_logger = Mock()
        logger = StructuredLogger('test', mock_logger)

        logger.info('info message')

        mock_logger.log.assert_called_once()
        call_args = mock_logger.log.call_args
        assert call_args[0][0] == logging.INFO

    def test_warning_logs_message(self):
        """Warning level logging."""
        from deployment.model_evaluation.logging_config import StructuredLogger

        mock_logger = Mock()
        logger = StructuredLogger('test', mock_logger)

        logger.warning('warning message')

        mock_logger.log.assert_called_once()
        call_args = mock_logger.log.call_args
        assert call_args[0][0] == logging.WARNING

    def test_error_logs_with_exc_info(self):
        """Error with exception info."""
        from deployment.model_evaluation.logging_config import StructuredLogger

        mock_logger = Mock()
        logger = StructuredLogger('test', mock_logger)

        logger.error('error message', exc_info=True)

        mock_logger.log.assert_called_once()
        call_args = mock_logger.log.call_args
        assert call_args[1].get('exc_info') is True

    def test_critical_logs_message(self):
        """Critical level logging."""
        from deployment.model_evaluation.logging_config import StructuredLogger

        mock_logger = Mock()
        logger = StructuredLogger('test', mock_logger)

        logger.critical('critical message')

        mock_logger.log.assert_called_once()
        call_args = mock_logger.log.call_args
        assert call_args[0][0] == logging.CRITICAL

    def test_with_extra_data(self):
        """Pass extra data with logging."""
        mock_logger = Mock()
        logger = lc.StructuredLogger('test', mock_logger)

        logger.info('message', extra={'key': 'value'})

        mock_logger.log.assert_called_once()
        call_args = mock_logger.log.call_args
        assert 'extra' in call_args[1]


class TestSetupLogging:
    """Tests for setup_logging function."""

    def test_setup_console_only(self):
        """Setup logging to console only."""
        logger = lc.setup_logging(console_level=logging.INFO)

        assert logger.level == logging.DEBUG
        assert len(logger.handlers) >= 1

    def test_setup_with_file(self):
        """Setup logging to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "test.log"

            logger = lc.setup_logging(log_file=log_file, console_level=logging.WARNING)

            assert log_file.exists()

    def test_setup_json_format(self):
        """Use JSON formatter for file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "test.json"

            logger = lc.setup_logging(
                log_file=log_file,
                console_level=logging.WARNING,
                log_format='json'
            )

            assert log_file.exists()

    def test_setup_text_format(self):
        """Use text formatter for file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "test.txt"

            logger = lc.setup_logging(
                log_file=log_file,
                console_level=logging.WARNING,
                log_format='text'
            )

            assert log_file.exists()

    def test_creates_log_directory(self):
        """Create log directory if not exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_dir = Path(tmpdir) / "nested" / "logs"
            log_file = log_dir / "test.log"

            logger = lc.setup_logging(log_file=log_file)

            assert log_file.exists()

    def test_setup_default_log_file(self):
        """Use default log file when log_dir provided."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_dir = Path(tmpdir)

            logger = lc.setup_logging(log_dir=log_dir)

            expected_file = log_dir / "evaluate.log"
            assert expected_file.exists()

    def test_multiple_handlers_cleared(self):
        """Clear existing handlers before adding new ones."""
        logger1 = lc.setup_logging(console_level=logging.INFO)
        initial_handlers = len(logger1.handlers)

        logger2 = lc.setup_logging(console_level=logging.WARNING)

        assert len(logger2.handlers) == initial_handlers


class TestLogContext:
    """Tests for LogContext class."""

    def test_temporary_level_change(self):
        """Temporarily change log level."""
        from deployment.model_evaluation.logging_config import LogContext

        test_logger = logging.getLogger('test_context')
        original_level = test_logger.level

        with LogContext(test_logger, logging.DEBUG):
            assert test_logger.level == logging.DEBUG

        assert test_logger.level == original_level

    def test_restores_original_level(self):
        """Restore original level on exit."""
        from deployment.model_evaluation.logging_config import LogContext

        test_logger = logging.getLogger('test_restore')
        test_logger.setLevel(logging.WARNING)

        with LogContext(test_logger, logging.DEBUG):
            pass

        assert test_logger.level == logging.WARNING


class TestLogFunctionCallDecorator:
    """Tests for log_function_call decorator."""

    def test_decorator_logs_call(self):
        """Log function call with args."""
        from deployment.model_evaluation.logging_config import log_function_call

        mock_logger = Mock()

        @log_function_call(mock_logger)
        def my_function(arg1, arg2):
            return arg1 + arg2

        my_function(1, 2)

        mock_logger.debug.assert_called()
        call_args = mock_logger.debug.call_args
        assert 'my_function' in str(call_args)

    def test_decorator_logs_result(self):
        """Log successful completion."""
        from deployment.model_evaluation.logging_config import log_function_call

        mock_logger = Mock()

        @log_function_call(mock_logger)
        def my_function(arg1, arg2):
            return arg1 + arg2

        my_function(1, 2)

        calls = mock_logger.debug.call_args_list
        assert any('success' in str(call) for call in calls)

    def test_decorator_logs_error(self):
        """Log error on exception."""
        from deployment.model_evaluation.logging_config import log_function_call

        mock_logger = Mock()

        @log_function_call(mock_logger)
        def failing_function():
            raise ValueError("test error")

        with pytest.raises(ValueError):
            failing_function()

        mock_logger.error.assert_called()
        call_args = mock_logger.error.call_args
        assert 'error' in str(call_args).lower()


class TestGetLogger:
    """Tests for get_logger function."""

    def test_get_logger_returns_structured_logger(self):
        """Return StructuredLogger instance."""
        logger = lc.get_logger('test_module')

        assert isinstance(logger, lc.StructuredLogger)
        assert logger.name == 'test_module'
